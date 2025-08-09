function [LHP] = price_LHP_t_student_optimized(discounts, recovery, Ku_list, K_d, corr, p, nu)
% ---------------------------------------------------------------------------------------------
% price_LHP_t_student_optimized - VERSIONE OTTIMIZZATA VETTORIALE
%
% INPUTS:
%   discounts – discount factor applied to the final payoff
%   recovery  – recovery rate on defaulted assets (e.g., 40%)
%   Ku_list   – vector of upper attachment points (e.g., [0.03, 0.06, 0.09, 0.12, 0.22])
%   K_d       – lower attachment point (typically 0 for equity tranche)
%   corr      – copula correlation parameter (ρ)
%   p         – marginal probability of default
%   nu        – degrees of freedom of the t-distribution
%
% OUTPUT:
%   LHP       – vector of discounted prices for all tranches under t-copula LHP model
% ---------------------------------------------------------------------------------------------

    % Step 1 – Compute the critical value in the t-distribution space
    K = tinv(p, nu);

    % Step 2 – Normalize tranche bounds in the cumulative loss space
    % Vectorize the upper bounds
    u_vec = Ku_list / (1 - recovery);  % vector of upper bounds
    d = K_d / (1 - recovery);          % scalar lower bound

    % Step 3 – Pre-compute constants
    sqrt_1_minus_corr = sqrt(1 - corr);
    sqrt_corr = sqrt(corr);
    
    % Step 4 – Define the vectorized integrand function
    % This function computes the contribution for ALL tranches simultaneously
    vectorized_integrand = @(z) compute_all_tranches_contribution(z, u_vec, d, K_d, Ku_list, ...
                                                                 sqrt_1_minus_corr, sqrt_corr, K, nu);

    % Step 5 – Single integration that computes expected losses for all tranches
    % The integrand returns a vector, so the integral returns a vector
    expected_losses = quadgk(vectorized_integrand, 0, 1, 'RelTol', 1e-10, 'AbsTol', 1e-12, ...
                            'ArrayValued', true);

    % Step 6 – Final prices: notional minus expected loss, discounted
    LHP = discounts * (1 - expected_losses);
end

function contributions = compute_all_tranches_contribution(z, u_vec, d, K_d, Ku_list, ...
                                                         sqrt_1_minus_corr, sqrt_corr, K, nu)
% ---------------------------------------------------------------------------------------------
% compute_all_tranches_contribution
%
% Computes the integrand contribution for all tranches simultaneously at a given z value.
% This eliminates the need for separate integrations for each tranche.
% ---------------------------------------------------------------------------------------------

    % Handle edge cases
    if z <= 1e-15 || z >= (1 - 1e-15)
        contributions = zeros(size(u_vec));
        return;
    end
    
    % Step 1 – Inverse mapping from cumulative loss z to systemic factor y
    % Handles the transformation y*(z) = (sqrt(1-ρ) * T⁻¹(z) - K) / sqrt(ρ)
    try
        tinv_z = tinv(z, nu);
        y_star = (sqrt_1_minus_corr * tinv_z - K) / sqrt_corr;
    catch
        contributions = zeros(size(u_vec));
        return;
    end
    
    % Step 2 – Compute transformed density (Jacobian of transformation)
    try
        tpdf_y_neg = tpdf(-y_star, nu);
        tpdf_tinv_z = tpdf(tinv_z, nu);
        
        if tpdf_tinv_z < 1e-15
            prob_z = 0;
        else
            prob_z = (sqrt_1_minus_corr / sqrt_corr) * tpdf_y_neg / tpdf_tinv_z;
        end
    catch
        contributions = zeros(size(u_vec));
        return;
    end
    
    % Step 3 – Vectorized tranche loss computation for all tranches
    % l(z) = min(max(((1-recovery)*z - K_d) / (K_u - K_d), 0), 1)
    numerator = (1 - recovery) * z - K_d;  % Scalar
    denominators = Ku_list - K_d;           % Vector
    
    % Vectorized tranche losses
    tranche_losses = min(max(numerator ./ denominators, 0), 1);
    
    % Step 4 – Final contributions for all tranches
    contributions = tranche_losses * prob_z;  % Element-wise multiplication with scalar
end