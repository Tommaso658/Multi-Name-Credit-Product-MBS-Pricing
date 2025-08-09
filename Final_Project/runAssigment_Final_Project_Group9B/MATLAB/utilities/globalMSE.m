function mse = globalMSE(nu, rho, discount, recovery, p, ...
                         Ku_list, Kd, rho_mkt)
% -------------------------------------------------------------------------
% globalMSE – Compute total Mean Squared Error between model-implied and 
% market-implied tranche prices using the LHP assumption for the t-student 
% copula model.
%
% INPUTS:
%   nu        – degrees of freedom of the t-Student copula
%   rho       – correlation parameter to be calibrated (same for all tranches)
%   discount  – discount factors (vector of zero-coupon bond prices)
%   recovery  – recovery rate assumed on defaulted exposures
%   p         – default probability over the horizon (e.g., 4 years)
%   Ku_list   – vector of upper attachment points for each tranche
%   Kd        – lower attachment point (common for all tranches in this setup)
%   rho_mkt   – vector of market-implied correlations for each tranche
%
% OUTPUT:
%   mse       – total mean squared error between model and market tranche prices
% -------------------------------------------------------------------------

    nTr = numel(Ku_list);  % Number of tranches to evaluate (one per Ku)
    
    mse = 0.0;              % Initialize the MSE accumulator

    % Loop over each tranche
    for j = 1:nTr
        Ku = Ku_list(j);    % Current upper attachment point (defines the tranche)
        
        % Compute the model-implied price of the tranche using the
        % t-student copula under the LHP approximation
        price_model = price_LHP_t_student(discount, recovery, Ku, Kd, ...
                                          rho, p, nu);
        
        % Compute the "market" price of the tranche using the Vasicek model
        % and the market-implied correlation for that specific tranche
        price_mkt   = price_LHP_vasicek(discount, recovery, Ku, Kd, ...
                                        rho_mkt(j), p);
        
        % Accumulate the squared error between model and market price
        mse = mse + (price_model - price_mkt)^2;
    end
end
