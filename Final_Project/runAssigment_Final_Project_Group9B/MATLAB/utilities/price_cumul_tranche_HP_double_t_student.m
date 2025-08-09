function [price] = price_cumul_tranche_HP_double_t_student(Kd, Ku, avg_recovery, rho, p, discount, I, nu)
% -----------------------------------------------------------------------------------------------------
% price_cumul_tranche_HP_double_t_student
%
% Computes the price of a cumulative tranche [Kd, Ku] under the double t-Student copula model
% for a finite number of exposures I, without using the LHP approximation.
%
% INPUTS:
%   Kd            – lower attachment point of the tranche (e.g., 0.03 for 3%)
%   Ku            – upper attachment point of the tranche
%   avg_recovery  – recovery rate on defaulted instruments (e.g., 0.4)
%   rho           – copula correlation parameter
%   p             – marginal default probability over the time horizon (e.g., 6%)
%   discount      – discount factor for present value computation
%   I             – number of exposures (finite portfolio, e.g., 500)
%   nu            – degrees of freedom for the t-Student distribution
%
% OUTPUT:
%   price         – tranche price (as a percentage of notional)
% -----------------------------------------------------------------------------------------------------

    % Step 1 – Normalize the tranche bounds in loss space [0, 1]
    d = Kd / (1 - avg_recovery);  % Lower threshold in normalized loss
    u = Ku / (1 - avg_recovery);  % Upper threshold in normalized loss

    % Step 2 – Define the tranche loss profile function
    % Maps number of defaults m to a loss proportion in the tranche,
    % clipped between 0 and 1, with linear increase between d and u
    loss_tranche = @(m) min(max((m / I - d) / (u - d), 0), 1);

    % Step 3 – Conditional default probability p(y)
    % Given systemic factor y ~ t(ν), this is the conditional default prob.
    p_y = @(y) tcdf((tinv(p, nu) - sqrt(rho) .* y) ./ sqrt(1 - rho), nu);

    % Step 4 – Probability of m defaults given systemic factor y
    % This follows a binomial distribution with parameter p_y(y)
    p_m_given_y = @(m, y) nchoosek(I, m) .* p_y(y).^m .* (1 - p_y(y)).^(I - m);

    % Step 5 – Compute expected loss by summing over m = 0 to I
    expected_loss = 0;  % Initialize expected loss accumulator

    for m = 0:I
        % For each number of defaults m:
        % Integrate over systemic factor y (with t-distribution PDF)
        expected_loss = expected_loss + ...
            quadgk(@(y) p_m_given_y(m, y) .* tpdf(y, nu), -6, 6) ...
            * loss_tranche(m);
    end

    % Step 6 – Compute the final tranche price:
    % price = discounted value of notional minus expected loss
    price = discount * (1 - expected_loss);
end
