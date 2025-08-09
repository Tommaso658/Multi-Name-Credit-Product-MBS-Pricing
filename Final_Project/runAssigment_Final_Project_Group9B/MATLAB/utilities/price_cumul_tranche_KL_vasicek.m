function price = price_cumul_tranche_KL_vasicek(kd, ku, recovery, rho, p, disc, I)
% -----------------------------------------------------------------------------------------
% price_cumul_tranche_KL_vasicek
%
% Computes the price of a cumulative tranche [kd, ku] under the Vasicek model
% using a Kullback-Leibler (KL) entropy-based analytical approximation.
%
% INPUTS:
%   kd       – lower attachment point of the tranche (in percentage of notional)
%   ku       – upper attachment point of the tranche
%   recovery – recovery rate on defaulted assets
%   rho      – asset correlation (Gaussian copula)
%   p        – default probability (at horizon)
%   disc     – discount factor
%   I        – number of exposures (finite, used in KL expansion)
%
% OUTPUT:
%   price    – approximated tranche price (percentage of notional)
% -----------------------------------------------------------------------------------------

    % Transform attachment points to loss space [0, 1]
    d = kd / (1 - recovery);  % normalized lower threshold
    u = ku / (1 - recovery);  % normalized upper threshold

    % Loss profile for the tranche: linear between d and u, bounded in [0, 1]
    loss_tr = @(z) max(0, min(1, (z - d) / (u - d)));

    % Conditional default probability given common factor y
    p_y = @(y) normcdf((norminv(p) - sqrt(rho) .* y) ./ sqrt(1 - rho));

    % Kullback-Leibler divergence between binomial and its conditional Bernoulli version
    K = @(z, y) z .* log(z ./ p_y(y)) + (1 - z) .* log((1 - z) ./ (1 - p_y(y)));

    % Normalization factor for saddlepoint approximation
    C1 = @(z) sqrt(I ./ (2 * pi .* z .* (1 - z)));

    % Integrand: joint density of (z, y) times the loss function and normal pdf
    integranda_zy = @(z, y) C1(z) .* exp(-I .* K(z, y)) .* loss_tr(z) .* normpdf(y);

    % Integration grid for y (common factor), from -6 to +6 standard deviations
    dy = 0.1;                      % discretization step for y
    discrete_y = (-6:dy:6)';       % grid points

    % Initialize integrand vector for outer integral (over y)
    integranda_y = zeros(length(discrete_y) - 1, 1);

    % Loop over y intervals and integrate over z ∈ [0,1] for each y midpoint
    for i = 1:(length(discrete_y) - 1)
        y_mid = (discrete_y(i) + discrete_y(i + 1)) / 2;
        integranda_y(i) = quadgk(@(z) integranda_zy(z, y_mid), 0, 1);
    end

    % Compute expected loss by integrating over y (weighted sum over y intervals)
    exp_loss = sum((discrete_y(2:end) - discrete_y(1:end-1)) .* integranda_y);

    % Correction: subtract the "atomic" term at full survival (no defaults)
    atom_fun = @(y) (1 - p_y(y)).^I .* normpdf(y);  % survival probability
    atom_term = ku * integral(atom_fun, -Inf, Inf, 'RelTol', 1e-8);  % small probability of no losses
    exp_loss = exp_loss - atom_term;

    % Final price of the tranche: notional minus expected loss, discounted
    price = disc * (1 - exp_loss);
end
