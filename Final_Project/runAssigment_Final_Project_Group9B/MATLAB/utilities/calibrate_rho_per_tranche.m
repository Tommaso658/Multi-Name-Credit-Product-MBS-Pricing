function [rho_matrix, MSE_nu, nu_star, rho_star_vec, nu_list_imp] = calibrate_rho_per_tranche(discounts, recovery, Ku_list, Kd_calibration, corr_mkt, p, nu_list_imp, optsFS)
%CALIBRATE_RHO_PER_TRANCHE  Calibrate tranche‑specific correlations under the
%double t‑Student copula (LHP approximation).
%
%   This routine iterates over a grid of degrees of freedom ν, and for each
%   tranche [Kd, Ku] solves an inverse pricing problem to find the asset
%   correlation ρ*(ν) that makes the model price of the tranche match the
%   Vasicek price obtained with the market‑implied correlation.  For every
%   ν the mean‑squared error (MSE) between calibrated correlations and
%   market correlations is computed; the function returns the ν* that
%   minimises this MSE together with the associated correlation vector.
%
%   INPUTS
%       discounts      – Discount‑factor curve (vector)
%       recovery       – Recovery rate (scalar)
%       Ku_list        – Row/column vector of upper attachment points (Ku)
%       Kd_calibration – Lower attachment point Kd used for every tranche
%       corr_mkt       – Vector of market‑implied tranche correlations ρ_mkt
%       p              – Default probability over the time horizon T
%       nu_list_imp    – (optional) Grid of ν values. Default linspace(2,151,151)
%       optsFS         – (optional) optimoptions for FSOLVE. Default: silent,
%                        TolFun=1e‑12, MaxIter=500
%
%   OUTPUTS
%       rho_matrix  – n_nu × n_tranches matrix with calibrated ρ*(ν) values
%       MSE_nu      – Vector of length n_nu with MSE(ν)
%       nu_star     – ν* that minimises the MSE
%       rho_star_vec– Row vector of calibrated ρ* at ν*
%       nu_list_imp – The ν grid actually used (useful when omitted)
%
%   EXAMPLE
%       nu_grid = linspace(2,151,151);
%       [rho_mat, MSE, nu_star, rho_star] = calibrate_rho_per_tranche( ...
%            discounts, 0.4, Ku_list, 0, corr_mkt, 0.06, nu_grid);
%
%   Author:  Group 9b – Financial Engineering Project AY 2024/2025
% %   Date  :  30‑May‑2025
%
%--------------------------------------------------------------------------

    % --------------------  Argument defaults  ---------------------------
    if nargin < 7 || isempty(nu_list_imp)
        nu_list_imp = linspace(2,151,151);  % Same grid as in the script
    end
    if nargin < 8 || isempty(optsFS)
        optsFS = optimoptions('fsolve', 'Display','off',               ...
                                        'TolFun',1e-12,              ...
                                        'MaxIterations',500);
    end

    % --------------------  Pre‑allocation  -----------------------------
    n_nu       = numel(nu_list_imp);
    n_tranches = numel(Ku_list);

    rho_matrix = NaN(n_nu, n_tranches);
    MSE_nu     = NaN(n_nu, 1);

    % --------------------  Grid search on ν  ---------------------------
    for i = 1:n_nu
        nu = nu_list_imp(i);

        for j = 1:n_tranches
            Ku = Ku_list(j);
            Kd = Kd_calibration;

            % Target price via Vasicek using market‑implied correlation
            price_target = price_LHP_vasicek(discounts, recovery, Ku, Kd,   ...
                                             corr_mkt(j), p);

            % Anonymous function: difference between t‑Student and target
            eq_fct = @(rho) price_LHP_t_student(discounts, recovery, Ku, Kd, ...
                                                rho, p, nu) - price_target;

            try
                rho_star = fsolve(eq_fct, 0.30, optsFS);  % Starting guess 0.30
            catch
                warning('calibrate\u005frho\u005fper\u005ftranche:NoConvergence', ...
                        'FSOLVE did not converge for ν = %g, tranche %d', nu, j);
                rho_star = NaN;  % Keep NaN so that MSE uses omitnan
            end

            rho_matrix(i, j) = rho_star;
        end

        % Compute MSE(ν) across tranches, safely ignoring NaN if any
        diff_sq    = (rho_matrix(i, :) - corr_mkt).^2;
        MSE_nu(i) = mean(diff_sq, 'omitnan');
    end

    % --------------------  Select best ν  ------------------------------
    [~, idx_best] = min(MSE_nu);
    nu_star       = nu_list_imp(idx_best);
    rho_star_vec  = rho_matrix(idx_best, :);
end
