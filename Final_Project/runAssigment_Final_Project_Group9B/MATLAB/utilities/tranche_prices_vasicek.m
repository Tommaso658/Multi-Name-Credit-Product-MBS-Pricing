function tranche_prices = tranche_prices_vasicek(I_list, Kd_calibration, Kd_list, Ku_list, recovery, corr_model, p, discount, flag)
% ------------------------------------------------------------------------------------------------------------
% tranche_prices_vasicek
%
% Computes the prices of multiple tranches [Kd, Ku] using the Vasicek copula model.
% The function allows three different pricing methods depending on the value of 'flag':
%   - flag = 0 : Exact pricing with finite number of exposures (binomial model)
%   - flag = 1 : KL approximation (analytical approach)
%   - flag = 2 : LHP approximation (Large Homogeneous Portfolio, I → ∞)
%
% INPUTS:
%   I_list          – vector of portfolio sizes (number of exposures)
%   Kd_calibration  – fixed lower attachment point for all tranches
%   Kd_list         – list of lower attachment points (not used here, may be redundant)
%   Ku_list         – list of upper attachment points (defining tranches)
%   recovery        – recovery rate (e.g., 0.4)
%   corr_model      – correlation parameter ρ for the Vasicek model
%   p               – marginal default probability
%   discount        – discount factor (applied uniformly)
%   flag            – pricing method selector:
%                       0 = Exact,
%                       1 = KL Approximation,
%                       2 = LHP Approximation
%
% OUTPUT:
%   tranche_prices  – matrix [length(I_list) × length(Ku_list)] of tranche prices
% ------------------------------------------------------------------------------------------------------------

    % ================
    % Case 1: Exact Pricing (finite I)
    % ================
    if flag == 0
        cumulative_prices_HP = zeros(length(I_list), length(Ku_list));  % Initialize result matrix

        for j = 1:length(I_list)  % Loop over different portfolio sizes
            for i = 1:length(Ku_list)  % Loop over tranches
                cumulative_prices_HP(j, i) = ...
                    price_cumul_tranche_HP_vasicek( ...
                        Kd_calibration, Ku_list(i), recovery, ...
                        corr_model, p, discount, round(I_list(j)));
            end
        end
        
        tranche_prices = cumulative_prices_HP;  % Assign to output

    end

    % ================
    % Case 2: KL Approximation (analytical)
    % ================
    if flag == 1
        cumulative_prices_KL = zeros(length(I_list), length(Ku_list));  % Initialize result matrix

        for j = 1:length(I_list)
            for i = 1:length(Ku_list)
                cumulative_prices_KL(j, i) = ...
                    price_cumul_tranche_KL_vasicek( ...
                        Kd_calibration, Ku_list(i), recovery, ...
                        corr_model, p, discount, round(I_list(j)));
            end
        end
        
        tranche_prices = cumulative_prices_KL;  % Assign to output
    end

    % ================
    % Case 3: LHP Approximation (I → ∞)
    % ================
    if flag == 2
        cumulative_prices_LHP = zeros(length(I_list), length(Ku_list));  % Initialize result matrix

        for j = 1:length(I_list)
            for i = 1:length(Ku_list)
                cumulative_prices_LHP(j, i) = ...
                    price_LHP_vasicek( ...
                        discount, recovery, Ku_list(i), ...
                        Kd_calibration, corr_model, p);
            end
        end
        
        tranche_prices = cumulative_prices_LHP;  % Assign to output
    end

end
