function [price] = price_cumul_tranche_HP_vasicek(Kd, Ku, avg_recovery, rho, p, discount, I)
% -----------------------------------------------------------------------------------------
% price_cumul_tranche_HP_vasicek
%
% Computes the price of a cumulative tranche [Kd, Ku] using the Vasicek model
% under a finite number of exposures I (i.e., no LHP approximation).
% This is a fully discrete model using conditional binomial distributions.
%
% INPUTS:
%   Kd            – lower attachment point of the tranche (e.g., 0.03 for 3%)
%   Ku            – upper attachment point of the tranche
%   avg_recovery  – average recovery rate on defaulted loans
%   rho           – correlation parameter in the Vasicek (Gaussian copula) model
%   p             – marginal probability of default (over time horizon)
%   discount      – discount factor applied to the payoff
%   I             – number of exposures (finite, no LHP approximation)
%
% OUTPUT:
%   price         – tranche price (% of notional), computed as 1 - expected loss
% -----------------------------------------------------------------------------------------

    % Step 1 – Normalize attachment points to loss space [0, 1]
    d = Kd / (1 - avg_recovery);  % loss threshold corresponding to Kd
    u = Ku / (1 - avg_recovery);  % loss threshold corresponding to Ku

    % Step 2 – Define tranche loss function
    % Given a number of defaults m out of I, maps to a percentage loss in the tranche
    % The function is piecewise linear between d and u, clipped to [0, 1]
    loss_tranche = @(m) min(max((m / I - d) / (u - d), 0), 1);

    % Step 3 – Define conditional default probability p(y)
    % For a given realization of the common factor y (standard normal),
    % the conditional probability that a name defaults is given by the probit form.
    p_y = @(y) normcdf((norminv(p) - sqrt(rho) .* y) ./ sqrt(1 - rho));

    % Step 4 – Define the conditional probability of m defaults given y
    % Follows a binomial distribution: P[m defaults | y] = Binomial(I, p_y(y))
    p_m_given_y = @(m, y) nchoosek(I, m) .* p_y(y).^m .* (1 - p_y(y)).^(I - m);

    % Step 5 – Initialize expected loss
    expected_loss = 0;

    % Step 6 – Loop over all possible default counts m = 0, ..., I
    % For each m, compute the expected value of the tranche loss given y, 
    % then integrate over y to get the unconditional expectation.
    for m = 0:I
        % Integrate over y using quadrature, assuming y ∈ [-6, 6]
        % normpdf(y) is the density of y ~ N(0,1)
        expected_loss = expected_loss + ...
            quadgk(@(y) p_m_given_y(m, y) .* normpdf(y), -6, 6) ...
            * loss_tranche(m);
    end

    % Step 7 – Compute final price: price = discounted expected payoff
    % The payoff is 1 - expected loss, since losses are subtracted from notional
    price = discount * (1 - expected_loss);
end
