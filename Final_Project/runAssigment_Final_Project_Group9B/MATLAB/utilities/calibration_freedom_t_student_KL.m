function [MSE_min, rho_model, idx_min, nu_list, MSE_list] = calibration_freedom_t_student_KL(discounts, corr_mkt, Kd_calibration, Ku_list, recovery, p, I, displayFlag)
% -------------------------------------------------------------------------------------------
% calibration_freedom_t_student_KL
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
    MSE_list = inf(1, length(nu_list));        % Initialize MSE for each ν to infinity
    MSE_min = inf;                             % Current best MSE
    idx_min = 1;                               % Index of ν with the best MSE
    rho_model = 0.3;                           % Initial guess for ρ (can be overwritten)

    % --- Compute the target price of the equity tranche using the Vasicek model
    %     This is the market reference value we want to match with our model
    price_equity_mkt = price_cumul_tranche_KL_vasicek(Kd_calibration, Ku_list(1), ...
                                                      recovery, corr_mkt(1), ...
                                                      p, discounts, I);

    % --- Setup solver options for fsolve
    %     High tolerance and max iterations, plus optional verbosity
    if nargin < 7
        displayFlag = 'off';                   % Default to silent mode if not specified
    end

    options = optimoptions('fsolve', ...
                           'TolFun',1e-12, ...
                           'Display',displayFlag, ...
                           'MaxIterations',500);

    % --- Loop through all ν values in the list
    for i = 1:length(nu_list)
        nu = nu_list(i);                       % Current ν being tested

        % --- Define the root-finding problem to match equity tranche price
        %     We look for a ρ that makes the t-student price = market price
        equation = @(rho) price_cumul_tranche_KL_double_t(Kd_calibration, Ku_list(1), ...
                                                          recovery, rho, p, ...
                                                          discounts, I, nu) ...
                          - price_equity_mkt;

        % --- Solve for ρ using fsolve
        try
            rho_star = fsolve(equation, 0.3, options);  % Initial guess: ρ = 0.3
        catch
            warning('Convergence problem for nu = %.2f', nu);
            continue;                                  % Skip this ν if solver fails
        end

        % --- Compute total MSE across the remaining tranches (excluding equity)
        mse = 0;
        for j = 2:length(Ku_list)
            % Price from double t-student copula (KL approximation)
            price_t = price_cumul_tranche_KL_double_t(Kd_calibration, Ku_list(j), ...
                                                      recovery, rho_star, p, ...
                                                      discounts, I, nu);

            % Reference market price using Vasicek copula
            price_v = price_cumul_tranche_KL_vasicek(Kd_calibration, Ku_list(j), ...
                                                     recovery, corr_mkt(j), p, ...
                                                     discounts, I);

            % Sum squared pricing errors
            mse = mse + (price_t - price_v)^2;
        end

        % --- Save the result for this ν
        MSE_list(i) = mse;

        % --- Update best result if current MSE is smaller than previous best
        if mse < MSE_min
            MSE_min = mse;          % Update minimum MSE
            idx_min = i;            % Save index of ν*
            rho_model = rho_star;   % Save best correlation found
        end
    end
end
