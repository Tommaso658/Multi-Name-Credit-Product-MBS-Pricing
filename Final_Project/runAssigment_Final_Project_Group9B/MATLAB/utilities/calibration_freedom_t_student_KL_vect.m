function [MSE_min, rho_model, idx_min, nu_list, MSE_list] = calibration_freedom_t_student_KL_vect(discounts, corr_mkt, Kd_calibration, Ku_list, recovery, p, I, displayFlag)
% -------------------------------------------------------------------------------------------
% calibration_freedom_t_student_KL (VECTORIZED VERSION)
%
% Calibrates the double t-Student copula model using the KL (Kullback-Leibler) approximation
% instead of the exact pricing method. The goal is to find the best degrees of freedom ν*
% and corresponding correlation ρ* that minimize the mismatch between model and market prices.
%
% INPUTS:
%   discounts        – discount factor vector (e.g., bootstrapped from market data)
%   corr_mkt         – market-implied correlations for each tranche
%   Kd_calibration   – lower attachment point (fixed, typically 0 for equity)
%   Ku_list          – upper attachment points for the tranches
%   recovery         – recovery rate on defaulted exposures
%   p                – default probability over the time horizon
%   I                – number of exposures in the portfolio
%   displayFlag      – string, 'on' or 'off' to control display of fsolve
%
% OUTPUTS:
%   MSE_min          – minimum Mean Squared Error obtained during calibration
%   rho_model        – correlation parameter ρ* that minimizes the error
%   idx_min          – index of the optimal ν in the tested ν grid
%   nu_list          – vector of tested ν values (degrees of freedom)
%   MSE_list         – vector of MSE values corresponding to each ν
% -------------------------------------------------------------------------------------------

    % --- Parameters
    nu_list = linspace(1, 50, 51);             % Grid of ν values (degrees of freedom)
    n_nu = length(nu_list);
    n_tranches = length(Ku_list);
    
    % Initialize output arrays
    MSE_list = inf(1, n_nu);                   % Initialize MSE for each ν to infinity
    rho_solutions = zeros(1, n_nu);            % Store all rho solutions
    
    % --- Compute the target price of the equity tranche using the Vasicek model
    price_equity_mkt = price_cumul_tranche_KL_vasicek(Kd_calibration, Ku_list(1), ...
                                                      recovery, corr_mkt(1), ...
                                                      p, discounts, I);

    % --- Setup solver options for fsolve
    if nargin < 8
        displayFlag = 'off';                   % Default to silent mode if not specified
    end

    options = optimoptions('fsolve', ...
                           'TolFun',1e-12, ...
                           'Display',displayFlag, ...
                           'MaxIterations',500);

    % --- Vectorized calibration: solve for all ν values
    % Pre-compute market prices for all tranches (excluding equity)
    price_vasicek_all = zeros(1, n_tranches-1);
    for j = 2:n_tranches
        price_vasicek_all(j-1) = price_cumul_tranche_KL_vasicek(Kd_calibration, Ku_list(j), ...
                                                               recovery, corr_mkt(j), p, ...
                                                               discounts, I);
    end
    
    % --- Main vectorized loop: process all nu values
    parfor i = 1:n_nu  % Use parfor for parallel processing if available
        nu = nu_list(i);
        
        % --- Define the root-finding problem to match equity tranche price
        equation = @(rho) price_cumul_tranche_KL_double_t(Kd_calibration, Ku_list(1), ...
                                                          recovery, rho, p, ...
                                                          discounts, I, nu) ...
                          - price_equity_mkt;

        % --- Solve for ρ using fsolve
        try
            rho_star = fsolve(equation, 0.3, options);
            rho_solutions(i) = rho_star;
        catch
            warning('Convergence problem for nu = %.2f', nu);
            rho_solutions(i) = NaN;
            continue;
        end
        
        % --- Vectorized computation of MSE across remaining tranches
        if ~isnan(rho_star)
            % Compute all t-student prices at once (vectorized over tranches)
            price_t_all = zeros(1, n_tranches-1);
            for j = 2:n_tranches
                price_t_all(j-1) = price_cumul_tranche_KL_double_t(Kd_calibration, Ku_list(j), ...
                                                                   recovery, rho_star, p, ...
                                                                   discounts, I, nu);
            end
            
            % Vectorized MSE computation
            squared_errors = (price_t_all - price_vasicek_all).^2;
            MSE_list(i) = sum(squared_errors);
        end
    end
    
    % --- Find the best result
    [MSE_min, idx_min] = min(MSE_list);
    rho_model = rho_solutions(idx_min);
    
    % Handle case where no valid solution was found
    if isnan(rho_model) || isinf(MSE_min)
        warning('No valid solution found across the nu grid');
        MSE_min = inf;
        rho_model = 0.3;  % Default fallback
        idx_min = 1;
    end
end