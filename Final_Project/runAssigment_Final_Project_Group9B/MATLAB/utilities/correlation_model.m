function [final_corr_model, min_idx, min_MSE] = correlation_model(nu_list, MSE_list, min_MSE, ...
                                                                   corr_mkt, Kd_calibration, ...
                                                                   Ku_list, recovery, p, discount)
% ------------------------------------------------------------------------------
% correlation_model – Find the optimal ν (degrees of freedom) that minimizes
% the MSE between model-implied and market-implied tranche correlations.
%
% INPUTS:
%   nu_list        – vector of ν values (degrees of freedom) to test
%   MSE_list       – vector of MSEs (will be updated inside the loop)
%   min_MSE        – current minimum MSE (starting value, may be updated)
%   corr_mkt       – vector of market-implied correlations per tranche
%   Kd_calibration – fixed lower attachment point for all tranches
%   Ku_list        – vector of upper attachment points (one per tranche)
%   recovery       – recovery rate assumed for each exposure
%   p              – default probability at maturity
%   discount       – discount factors (e.g., from zero curve)
%
% OUTPUTS:
%   final_corr_model – vector of model-implied correlations (best fit)
%   min_idx          – index in nu_list corresponding to min MSE
%   min_MSE          – value of the minimum MSE found
% ------------------------------------------------------------------------------

    min_idx = 1;  % Initialize index of minimum MSE (will be updated if a better value is found)

    % Loop through each value of ν in the candidate list
    for i = 1:length(MSE_list)

        % For the current ν, perform calibration:
        % - Find the model-implied correlation that best matches market prices
        % - Compute the corresponding MSE
        [MSE_list(i), corr_model] = calibration_freedom(discount, ...
                                                        corr_mkt, ...
                                                        Kd_calibration, ...
                                                        Ku_list, ...
                                                        recovery, ...
                                                        p, ...
                                                        nu_list(i));

        % Check if this MSE is better (lower) than the current minimum
        if MSE_list(i) < min_MSE
            min_MSE = MSE_list(i);        % Update minimum MSE
            min_idx = i;                  % Save index of optimal ν
            final_corr_model = corr_model;  % Save corresponding calibrated correlation(s)
        end
    end
end
S