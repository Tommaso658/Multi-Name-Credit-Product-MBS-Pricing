function [MSE_min, rho_model, idx_min, nu_grid, nu_star, rho_star_vec, MSE_list] = calibration_t_student( ...
            discount, corr_mkt, Kd_calibration, Ku_list, recovery, p, displayFlag, flag)

% -------------------------------------------------------------------------
% Purpose:
%   Calibrate the t-Student copula model to market tranche correlations
%
% Inputs:
%   discount         : Discount factor (scalar)
%   corr_mkt         : Market-implied asset correlations for each tranche (1 x N vector)
%   Kd_calibration   : Lower attachment point used for calibration (scalar)
%   Ku_list          : Upper detachment points of tranches (1 x N vector)
%   recovery         : Recovery rate (scalar)
%   p                : Default probability (scalar)
%   displayFlag      : Controls display options (currently unused here)
%   flag             : Calibration mode selector:
%                      0 = ν fixed, calibrate ρ (calls external function)
%                      1 = calibrate ρ for each ν, then select ν* (main focus here)
%                      2 = global search over (ν, ρ)
%
% Outputs:
%   MSE_min          : Minimum Mean Squared Error across all ν (or from global opt)
%   rho_model        : (Only for flag=0) Calibrated ρ
%   idx_min          : Index of ν* in the ν grid
%   nu_grid          : List of ν values tested (only for flag=1)
%   nu_star          : Optimal degrees of freedom ν*
%   rho_star_vec     : Vector of calibrated ρ*(ν*) for each tranche
%   MSE_list         : MSE associated with each ν in the grid

% Default outputs
MSE_min = NaN; rho_model = NaN; idx_min = NaN;
nu_star = NaN; rho_star_vec = NaN; MSE_list = NaN; nu_grid = NaN;

% -------------------------------------------------------------------------
if flag == 0
    % --- Simple calibration with fixed ν, optimize ρ for equity tranche only ---
    [MSE_min, rho_model, idx_min, nu_grid, MSE_list] = calibration_freedom_t_student_LHP( ...
        discount, corr_mkt, Kd_calibration, Ku_list, recovery, p, displayFlag);
    
% -------------------------------------------------------------------------
elseif flag == 1
    % --- Full calibration: calibrate ρ per tranches, then choose ν* minimizing MSE ---
    
    nu_grid = linspace(2, 151, 151);     % ν values from 2 to 151
    n_nu = numel(nu_grid);               % Number of ν values
    n_tranches = numel(Ku_list);         % Number of tranches

    rho_matrix = NaN(n_nu, n_tranches);  % Matrix to store calibrated ρ*(ν, tranche)
    MSE_list = NaN(n_nu, 1);             % MSE at each ν

    % fsolve options for solving ρ(ν) per tranche
    optsFS = optimoptions('fsolve', 'Display', 'off', ...
                                     'TolFun', 1e-12, ...
                                     'MaxIterations', 500);

    % Loop over ν values
    for i = 1:n_nu
        nu = nu_grid(i);

        for j = 1:n_tranches
            Ku = Ku_list(j);
            Kd = Kd_calibration;

            % Market-implied price under Vasicek model (used as target)
            price_target = price_LHP_vasicek(discount, recovery, Ku, Kd, ...
                                             corr_mkt(j), p);

            % Define zero function: price_model_t_student(ρ, ν) - price_target
            eq_fct = @(rho) price_LHP_t_student(discount, recovery, Ku, Kd, ...
                                                rho, p, nu) - price_target;

            try
                rho_star = fsolve(eq_fct, 0.30, optsFS);  % Initial guess: 0.30
            catch
                warning('fsolve did not converge for ν = %.1f, tranche %d', nu, j);
                rho_star = NaN;
            end

            rho_matrix(i, j) = rho_star;
        end

        % Compute MSE between model-calibrated ρ and market correlations
        MSE_list(i) = mean((rho_matrix(i,:) - corr_mkt).^2, 'omitnan');
    end

    % Find ν* minimizing the MSE
    [MSE_min, idx_min] = min(MSE_list);
    nu_star = nu_grid(idx_min);
    rho_star_vec = rho_matrix(idx_min, :);

% -------------------------------------------------------------------------
elseif flag == 2
    % --- Global 2D optimization: search over (ν, ρ) jointly ---
    
    obj = @(x) globalMSE(x(1), x(2), discount, recovery, p, ...
                         Ku_list, Kd_calibration, corr_mkt);  % x = [ν, ρ]

    x0 = [10, 0.30];                      % Initial guess: [ν, ρ]
    opts = optimset('Display', 'iter');  % Display optimization steps

    [x_star, MSE_star] = fminsearch(obj, x0, opts);

    nu_star = x_star(1);
    rho_model = x_star(2);
    MSE_min = MSE_star;

    fprintf('ν* = %.4f   ρ* = %.4f   MSE = %.4e\n', nu_star, rho_model, MSE_min);

end
end
