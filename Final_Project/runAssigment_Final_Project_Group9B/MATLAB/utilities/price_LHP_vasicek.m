function [LHP] = price_LHP_vasicek(discounts, recovery, K_u, K_d, corr, p)
% ----------------------------------------------------------------------------------------
% price_LHP_vasicek
%
% Computes the price of a tranche [K_d, K_u] using the Vasicek model under
% the Large Homogeneous Portfolio (LHP) assumption.
%
% INPUTS:
%   discounts – discount factor for present value
%   recovery  – recovery rate (e.g., 40%)
%   K_u       – upper attachment point of the tranche (e.g., 0.06)
%   K_d       – lower attachment point of the tranche (e.g., 0.03)
%   corr      – asset correlation (ρ)
%   p         – marginal probability of default (PD)
%
% OUTPUT:
%   LHP       – tranche price under the Vasicek LHP model
% ----------------------------------------------------------------------------------------

    % Compute the default threshold in standard normal space
    % It is the quantile such that P(X ≤ K) = p, with X ~ N(0,1)
    K = norminv(p);

    % Normalize tranche bounds in loss space (scale up due to LGD = 1 - recovery)
    u = K_u / (1 - recovery);  % upper limit in loss domain
    d = K_d / (1 - recovery);  % lower limit in loss domain

    % Compute the value of systemic factor y* that implies a conditional default probability = z
    % Derived from solving Φ((K - sqrt(1 - ρ) Φ⁻¹(z)) / sqrt(ρ)) = z
    y_star = @(z) ((sqrt(1 - corr)) .* norminv(z) - K) ./ sqrt(corr);

    % Compute the derivative dz/dy using the change-of-variable formula
    % This gives the density of z under the transformation from y
    prob = @(z) (sqrt(1 - corr) / sqrt(corr)) .* normpdf(-y_star(z)) ./ normpdf(norminv(z));

    % Define the loss in the tranche as a function of cumulative portfolio loss z
    % It’s clipped between 0 and 1 and scales linearly between K_d and K_u
    l = @(z) min(max(((1 - recovery) .* z - K_d) ./ (K_u - K_d), 0), 1);

    % Define the integrand function f(z) = loss(z) * density(z)
    f = @(z) l(z) .* prob(z);

    % Compute the expected loss over z ∈ [0,1] via numerical integration
    ELHP = quadgk(f, 0, 1);  % Expected tranche loss under LHP

    % Final price = discounted notional * (1 - expected loss)
    LHP = discounts * (1 - ELHP);
end
