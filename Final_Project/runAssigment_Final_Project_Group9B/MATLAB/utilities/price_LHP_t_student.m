function [LHP] = price_LHP_t_student(discounts, recovery, K_u, K_d, corr, p, nu)
% ---------------------------------------------------------------------------------------------
% price_LHP_t_student
%
% Computes the price of a tranche [K_d, K_u] using the double t-Student copula model
% under the Large Homogeneous Portfolio (LHP) approximation.
%
% INPUTS:
%   discounts – discount factor applied to the final payoff
%   recovery  – recovery rate on defaulted assets (e.g., 40%)
%   K_u       – upper attachment point of the tranche (e.g., 6%)
%   K_d       – lower attachment point of the tranche (e.g., 3%)
%   corr      – copula correlation parameter (ρ)
%   p         – marginal probability of default
%   nu        – degrees of freedom of the t-distribution
%
% OUTPUT:
%   LHP       – discounted price of the tranche under the t-copula LHP model
% ---------------------------------------------------------------------------------------------

    % Step 1 – Compute the critical value in the t-distribution space
    % The threshold such that P(T ≤ K) = p, with T ~ t_ν
    K = tinv(p, nu);

    % Step 2 – Normalize tranche bounds in the cumulative loss space
    u = K_u / (1 - recovery);  % upper bound in normalized loss space
    d = K_d / (1 - recovery);  % lower bound in normalized loss space

    % Step 3 – Inverse mapping from cumulative loss z to the systemic factor y
    % Solves the equation: tcdf((K - sqrt(1 - ρ) * T⁻¹(z)) / sqrt(ρ), ν) = z
    y_star = @(z) ((sqrt(1 - corr)) .* tinv(z, nu) - K) ./ sqrt(corr);

    % Step 4 – Compute transformed density (dz/dy), i.e. the Jacobian of the transformation
    % This is the derivative of the mapping from y to z, in terms of t-distributions
    prob = @(z) (sqrt(1 - corr) / sqrt(corr)) .* ...
                tpdf(-y_star(z), nu) ./ tpdf(tinv(z, nu), nu);

    % Step 5 – Define tranche loss as a function of z (cumulative portfolio loss)
    % Linear increase from 0 to 1 between d and u, clipped between [0, 1]
    l = @(z) min(max(((1 - recovery) .* z - K_d) ./ (K_u - K_d), 0), 1);

    % Step 6 – Define the integrand: tranche loss × density under the transformed variable
    f = @(z) l(z) .* prob(z);

    % Step 7 – Compute expected loss in the tranche via integration over z ∈ [0, 1]
    ELHP = quadgk(f, 0, 1);

    % Step 8 – Final price: notional minus expected loss, discounted
    LHP = discounts * (1 - ELHP);
end
