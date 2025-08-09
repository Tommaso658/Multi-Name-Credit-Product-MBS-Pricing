function tranche_prices = tranche_prices_double_t(nu_optimal, I_list, Kd_calibration, Kd_list, Ku_list, recovery, corr_model, p, discount, flag)
% -------------------------------------------------------------------------------------------------------
% tranche_prices_double_t
%
% Computes the prices of a set of tranches under the double t-Student copula model.
% It supports three modes:
%   - Exact pricing using simulation over finite I (flag = 0)
%   - KL approximation for analytical pricing (flag = 1)
%   - Large Homogeneous Portfolio (LHP) approximation (flag = 2)
%
% INPUTS:
%   nu_optimal      – degrees of freedom of the t-distribution (ν*)
%   I_list          – vector of exposure counts (e.g., [10, 20, ..., 1000])
%   Kd_calibration  – fixed lower attachment point for all tranches (e.g., 0)
%   Kd_list         – vector of K_d (may be unused in this version)
%   Ku_list         – vector of K_u, defining tranche upper bounds
%   recovery        – recovery rate on defaulted assets
%   corr_model      – calibrated copula correlation ρ*
%   p               – marginal default probability
%   discount        – discount factor to apply to final price
%   flag            – selector for method: 
%                      0 = exact pricing,
%                      1 = KL approximation,
%                      2 = LHP approximation
%
% OUTPUT:
%   tranche_prices  – matrix of tranche prices: rows = I values, columns = tranches
% -------------------------------------------------------------------------------------------------------

    % ================
    % Case 1: Real Price (Exact Simulation)
    % ================
    if flag == 0
        % Initialize the result matrix [rows = I values, columns = Ku tranches]
        cumulative_prices_HP = zeros(length(I_list), length(Ku_list));

        % For each value of I (portfolio size)
        for j = 1:length(I_list)
            for i = 1:length(Ku_list)
                % Compute the tranche price using finite I and t-copula model
                cumulative_prices_HP(j, i) = price_cumul_tranche_HP_double_t_student( ...
                    Kd_calibration, Ku_list(i), recovery, ...
                    corr_model, p, discount, round(I_list(j)), nu_optimal);
            end
        end

        tranche_prices = cumulative_prices_HP;  % Output matrix

    end

    % ================
    % Case 2: KL Approximation (Analytical)
    % ================
    if flag == 1
        cumulative_prices_KL = zeros(length(I_list), length(Ku_list));

        for j = 1:length(I_list)
            for i = 1:length(Ku_list)
                % Use KL-based pricing formula (approximate but fast)
                cumulative_prices_KL(j, i) = price_cumul_tranche_KL_double_t( ...
                    Kd_calibration, Ku_list(i), recovery, ...
                    corr_model, p, discount, round(I_list(j)), nu_optimal);
            end
        end

        tranche_prices = cumulative_prices_KL;  % Output matrix
    end

    % ================
    % Case 3: LHP Approximation (Limit case I → ∞)
    % ================
    if flag == 2
        cumulative_prices_LHP = zeros(length(I_list), length(Ku_list));

        for j = 1:length(I_list)
            for i = 1:length(Ku_list)
                % Uses analytical expression assuming an infinite pool (LHP)
                cumulative_prices_LHP(j, i) = price_LHP_t_student( ...
                    discount, recovery, Ku_list(i), ...
                    Kd_calibration, corr_model, p, nu_optimal);
            end
        end

        tranche_prices = cumulative_prices_LHP;  % Output matrix
    end
end
