function [MSE_min, rho_model, idx_min, nu_list, MSE_list] = calibration_freedom_t_student_LHP(discounts, corr_mkt, Kd_calibration, Ku_list, recovery, p, displayFlag)
% -----------------------------------------------------------------------------------------
% calibration_freedom_t_student_LHP
%
% Calibrates the double t-Student copula model under the Large Homogeneous Portfolio (LHP)
% assumption by finding the optimal degrees of freedom (ν) that minimize the pricing error.
%
% INPUTS:
%   discounts        – vector of discount factors (from bootstrap or given)
%   corr_mkt         – vector of market-implied correlations (Vasicek-based)
%   Kd_calibration   – fixed lower attachment point for tranches (typically 0 for equity)
%   Ku_list          – vector of upper attachment points for tranches
%   recovery         – recovery rate (e.g., 40%)
%   p                – default probability over time horizon (e.g., 6% over 4 years)
%   displayFlag      – 'on' or 'off' to control fsolve verbosity
%
% OUTPUTS:
%   MSE_min          – minimum mean squared error achieved
%   rho_model        – corresponding correlation parameter (ρ*) for best ν*
%   idx_min          – index in nu_list of optimal ν*
%   nu_list          – vector of tested ν values (degrees of freedom)
%   MSE_list         – vector of MSEs associated with each ν in nu_list
% -----------------------------------------------------------------------------------------

    % --- Parameters
    nu_list = linspace(1, 50, 50);            % Candidate degrees of freedom (from 1 to 50)
    MSE_list = inf(1, length(nu_list));       % Initialize all MSE values to ∞
    MSE_min = inf;                            % Minimum MSE so far
    idx_min = 1;                              % Index of ν* giving the minimum MSE
    rho_model = 0.3;                          % Initial value for ρ* (can be overwritten later)

    % --- Market price of the equity tranche computed using the Vasicek model
    %     Used as the target for calibration
    price_equity_mkt = price_LHP_vasicek(discounts, recovery, Ku_list(1), ...
                                         Kd_calibration, corr_mkt(1), p);

    % --- Setup fsolve options
    %     Use high precision and disable verbose output unless specified
    if nargin < 7
        displayFlag = 'off';                  % Default display mode
    end

    options = optimoptions('fsolve', ...
                           'TolFun',1e-12, ...    % High precision tolerance
                           'Display',displayFlag, ...
                           'MaxIterations',500);  % Ensure deep search in optimization

    % --- Loop over all candidate degrees of freedom (ν)
    for i = 1:length(nu_list)
        nu = nu_list(i);                      % Current degrees of freedom to test

        % --- Define the equation to solve: price_t(ρ) - price_equity_mkt = 0
        %     Goal: find ρ such that model price (t-copula) matches market price
        equation = @(rho) price_LHP_t_student(discounts, recovery, Ku_list(1), ...
                                              Kd_calibration, rho, p, nu) - price_equity_mkt;

        % --- Run the optimization to solve for ρ* at fixed ν
        try
            rho_star = fsolve(equation, 0.3, options);  % Initial guess for ρ = 0.3
        catch
            warning('Convergence problem for nu = %.2f', nu);  % Handle failed convergence
            continue;
        end

        % --- Compute total MSE over the REMAINING tranches (j = 2 to end)
        mse = 0;
        for j = 2:length(Ku_list)
            % Price using t-student model with calibrated ρ*
            price_t = price_LHP_t_student(discounts, recovery, Ku_list(j), ...
                                          Kd_calibration, rho_star, p, nu);

            % Price using Vasicek model and market-implied ρ for tranche j
            price_v = price_LHP_vasicek(discounts, recovery, Ku_list(j), ...
                                        Kd_calibration, corr_mkt(j), p);

            % Add squared pricing error to MSE
            mse = mse + (price_t - price_v)^2;
        end

        % --- Save result for current ν
        MSE_list(i) = mse;

        % --- Update if this is the best ν so far
        if mse < MSE_min
            MSE_min = mse;          % Best error so far
            idx_min = i;            % Best index (best ν*)
            rho_model = rho_star;   % Best correlation found
        end
    end
end
