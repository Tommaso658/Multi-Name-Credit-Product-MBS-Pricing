function tranche_prices = tranche_prices_vasicek_D(I_list, Kd_calibration, Kd_list, Ku_list, recovery, corr_model, p, discount, flag)
% ------------------------------------------------------------------------------------------------------------
% tranche_prices_vasicek
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
                price_cumul_tranche_HP_vasicek_optimized( ...
                    Kd_calibration, Ku_list(i), recovery, ...
                    corr_model, p, discount, round(I_list(j)));
        end
    end
    
    tranche_prices = cumulative_prices_HP;  % Assign to output

end

end

