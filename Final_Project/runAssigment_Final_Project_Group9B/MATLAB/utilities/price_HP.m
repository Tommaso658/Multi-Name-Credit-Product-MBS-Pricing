function [price_LHP, ELHP] = price_HP(discounts, dates, new_dates, I, recovery, K_u, K_d, corr, p)
% --------------------------------------------------------------------------------------------------
% price_HP
%
% Computes the price and expected loss of a tranche [K_d, K_u] under a homogeneous portfolio model
% using a KL-based analytical approximation (saddlepoint method).
%
% INPUTS:
%   discounts    – vector of discount factors corresponding to 'dates'
%   dates        – vector of original dates
%   new_dates    – date(s) to which the pricing refers (e.g. maturity date)
%   I            – number of exposures in the portfolio
%   recovery     – recovery rate on defaults
%   K_u          – upper attachment point of the tranche
%   K_d          – lower attachment point of the tranche
%   corr         – asset correlation (ρ)
%   p            – default probability (scalar)
%
% OUTPUTS:
%   price_LHP    – discounted tranche price (percentage of notional)
%   ELHP         – expected loss in the tranche (fraction of notional)
% --------------------------------------------------------------------------------------------------

    % Step 1 – Compute the default threshold in standard normal space
    K = norminv(p);  % Threshold such that P(X < K) = p, with X ~ N(0,1)

    % Step 2 – Extract the discount factor at maturity
    % Uses getDiscount() to align dates and extract the final value
    end_discounts = getDiscount(discounts, new_dates, dates, new_dates(1));
    end_discount = end_discounts(end);  % use the last value as maturity discount

    % Step 3 – Define the tranche loss function
    % Converts cumulative loss z into a proportion of loss in the tranche
    l = @(z) min(max(((1 - recovery) .* z - K_d) ./ (K_u - K_d), 0), 1);

    % Step 4 – Saddlepoint normalization factor (from large deviations theory)
    C = @(z) sqrt(I ./ (2 .* pi .* z .* (1 - z)));

    % Step 5 – Conditional default probability (simplified constant form)
    % May be meant as an approximation; this is not a function of y here
    p = @(y) normcdf((K - sqrt(corr)) ./ sqrt(1 - corr), 0, 1);

    % Step 6 – Kullback-Leibler divergence between z and p(y)
    % Measures the divergence between actual loss z and expected default p(y)
    K_bis = @(z, y) z .* log(z ./ p(y)) + ...
                    (1 - z) .* log((1 - z) ./ (1 - p(y)));

    % Step 7 – Full joint integrand over z and y
    % Includes: KL term, density over y, saddlepoint term, and tranche loss
    f = @(z, y) normpdf(y, 0, 1) .* C(z) .* ...
                exp(-I .* K_bis(z, p(y))) .* l(z);

    % Step 8 – Outer integral over y ∈ [−6, 6] of the inner integral in z ∈ [0, 1]
    inner_integral = @(y) quadgk(@(z) f(z, y), 0, 1);

    % Step 9 – Compute expected loss by integrating over y (common factor)
    ELHP = quadgk(inner_integral, -6, 6);

    % Step 10 – Compute final tranche price: notional minus expected loss, discounted
    price_LHP = end_discount * (1 - ELHP);
end
