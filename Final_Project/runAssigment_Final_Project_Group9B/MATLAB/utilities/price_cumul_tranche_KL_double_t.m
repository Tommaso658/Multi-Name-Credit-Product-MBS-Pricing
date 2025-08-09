function price = price_cumul_tranche_KL_double_t(kd, ku, recovery, rho, p, discount, I, nu)
% -----------------------------------------------------------------------------------------
% price_cumul_tranche_KL_double_t
%
% Computes the price of a tranche [kd, ku] using a double t-Student copula model
% with a Kullback-Leibler entropy-based approximation.
% The approach avoids Monte Carlo simulation and uses an analytical saddlepoint method.
%
% INPUTS:
%   kd        – lower attachment point of the tranche (e.g., 0.03 for 3%)
%   ku        – upper attachment point of the tranche
%   recovery  – recovery rate assumed on defaulted exposures
%   rho       – asset correlation (copula parameter)
%   p         – marginal default probability over the time horizon (e.g., 6%)
%   discount  – discount factor to present value the payoff
%   I         – number of exposures in the portfolio
%   nu        – degrees of freedom for the t-Student copula
%
% OUTPUT:
%   price     – approximated tranche price (% of notional)
% -----------------------------------------------------------------------------------------

    % Normalize the attachment points with respect to total loss (1 - recovery)
    d = kd / (1 - recovery);  % lower threshold in loss space
    u = ku / (1 - recovery);  % upper threshold in loss space

    % Tranche loss profile function:
    % Linearly increases from 0 to 1 between d and u, clipped to [0, 1]
    loss_tr = @(z) max(0, min(1, (z - d) / (u - d)));

    % Conditional default probability given common factor y (t-distributed)
    % This reflects the dependency structure induced by the t-copula
    p_y = @(y) tcdf((tinv(p, nu) - sqrt(rho) .* y) ./ sqrt(1 - rho), nu);

    % Kullback-Leibler divergence between Bernoulli(z) and conditional Bernoulli(p_y)
    % Measures the "distance" in entropy between the model and the target
    K = @(z, y) z .* log(z ./ p_y(y)) + (1 - z) .* log((1 - z) ./ (1 - p_y(y)));

    % Normalization factor in saddlepoint approximation
    C1 = @(z) sqrt(I ./ (2 * pi .* z .* (1 - z)));

    % Full integrand: joint density for (z, y) times the loss function and t pdf
    integranda_zy = @(z, y) C1(z) .* exp(-I .* K(z, y)) .* loss_tr(z) .* tpdf(y, nu);

    % Grid for y (the systemic risk factor), using midpoints between -6 and 6
    dy = 0.1;                          
    discrete_y = (-6:dy:6)';           

    % Initialize integral result per y slice
    integranda_y = zeros(length(discrete_y) - 1, 1);

    % Loop over y intervals and integrate over z ∈ [0, 1] at each y midpoint
    for i = 1:(length(discrete_y) - 1)
        y_mid = (discrete_y(i) + discrete_y(i + 1)) / 2;  % midpoint of y interval
        integranda_y(i) = quadgk(@(z) integranda_zy(z, y_mid), 0, 1);  % integrate over z
    end

    % Compute expected loss: integrate over y using trapezoidal rule approximation
    exp_loss = sum((discrete_y(2:end) - discrete_y(1:end-1)) .* integranda_y);

    % Atomic correction term: accounts for the probability of full survival (no defaults)
    % Underestimation can occur near the boundaries (especially for equity)
    atom_fun = @(y) (1 - p_y(y)).^I .* tpdf(y, nu);  % full survival probability * density
    atom_term = ku * quadgk(@(y) atom_fun(y), -30, 30);  % integrate over y domain
    exp_loss = exp_loss - atom_term;  % subtract from expected loss

    % Final tranche price: 1 - expected loss, discounted
    price = discount * (1 - exp_loss);
end
