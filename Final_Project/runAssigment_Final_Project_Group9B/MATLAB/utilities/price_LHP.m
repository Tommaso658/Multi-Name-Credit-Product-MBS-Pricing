function LHP = price_LHP(discount, recovery, K_u, K_d, corr, p)
% -------------------------------------------------------------------------------------
% price_LHP
%
% Approximates the price of a tranche [K_d, K_u] under a simplified
% Vasicek-style Large Homogeneous Portfolio (LHP) model.
% Note: This version differs slightly from standard Vasicek formulas due to the use
% of (1 - sqrt(corr)) instead of sqrt(1 - corr) – may reflect an approximation.
%
% INPUTS:
%   discount  – discount factor applied to expected cash flows
%   recovery  – recovery rate (e.g., 0.4)
%   K_u       – upper attachment point of the tranche
%   K_d       – lower attachment point of the tranche
%   corr      – asset correlation parameter (ρ)
%   p         – marginal default probability
%
% OUTPUT:
%   LHP       – tranche price (% of notional), under LHP approximation
% -------------------------------------------------------------------------------------

    % Step 1 – Compute the threshold K in normal space (quantile of default)
    K = norminv(p);  % critical value such that P(X ≤ K) = p, X ~ N(0,1)

    % Step 2 – Define inverse transformation y*(z) from portfolio loss z to systemic factor
    % Note: This form uses (1 - sqrt(ρ)) instead of standard sqrt(1 - ρ), which changes scaling
    y_star = @(z) -((1 - sqrt(corr)) .* norminv(z) - K) ./ sqrt(corr);

    % Step 3 – Compute the density adjustment (Jacobian of change of variables)
    % Derivative dz/dy, adapted to the transformation used above
    prob = @(z) ((1 - sqrt(corr)) / sqrt(corr)) .* ...
                normpdf(-y_star(z)) ./ normpdf(norminv(z));

    % Step 4 – Define the tranche loss function over cumulative loss z ∈ [0,1]
    % Linear increase between K_d and K_u, clipped to [0,1]
    l = @(z) min(max(((1 - recovery) .* z - K_d) ./ (K_u - K_d), 0), 1);

    % Step 5 – Integrand: loss function weighted by transformed density
    f = @(z) l(z) .* prob(z);

    % Step 6 – Compute expected loss by integrating f(z) over z ∈ [0,1]
    ELHP = quadgk(f, 0, 1);

    % Step 7 – Compute final price = discounted notional × (1 − expected loss)
    LHP = discount .* (1 - ELHP);
end
