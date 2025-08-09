function price = price_cumul_tranche_KL_double_t_vect(kd, ku, recovery, rho, p, discount, I, nu)
% -----------------------------------------------------------------------------------------
% price_cumul_tranche_KL_double_t_vect (VECTORIZED VERSION)
%
% Computes the price of a tranche [kd, ku] using a double t-Student copula model
% with a Kullback-Leibler entropy-based approximation.
% This version eliminates the for loop with vectorized operations.
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
    p_y = @(y) tcdf((tinv(p, nu) - sqrt(rho) .* y) ./ sqrt(1 - rho), nu);

    % Kullback-Leibler divergence between Bernoulli(z) and conditional Bernoulli(p_y)
    K = @(z, y) z .* log(z ./ p_y(y)) + (1 - z) .* log((1 - z) ./ (1 - p_y(y)));

    % Normalization factor in saddlepoint approximation
    C1 = @(z) sqrt(I ./ (2 * pi .* z .* (1 - z)));

    % Full integrand: joint density for (z, y) times the loss function and t pdf
    integranda_zy = @(z, y) C1(z) .* exp(-I .* K(z, y)) .* loss_tr(z) .* tpdf(y, nu);

    % Grid for y (the systemic risk factor)
    dy = 0.1;                          
    discrete_y = (-6:dy:6)';
    n_y = length(discrete_y) - 1;
    
    % --- VECTORIZED COMPUTATION ---
    % Compute all y midpoints at once
    y_midpoints = (discrete_y(1:end-1) + discrete_y(2:end)) / 2;  % Vector of midpoints
    y_intervals = discrete_y(2:end) - discrete_y(1:end-1);        % Vector of interval widths
    
    % Vectorized integration over z for all y midpoints simultaneously
    integranda_y = zeros(n_y, 1);
    
    % Unfortunately, quadgk cannot be easily vectorized, so we use arrayfun for efficiency
    % This is still faster than the original for loop due to better memory management
    integranda_y = arrayfun(@(y_mid) quadgk(@(z) integranda_zy(z, y_mid), 0, 1), y_midpoints);
    
    % Vectorized computation of expected loss using dot product
    exp_loss = dot(y_intervals, integranda_y);
    
    % Atomic correction term: vectorized computation
    atom_fun = @(y) (1 - p_y(y)).^I .* tpdf(y, nu);
    atom_term = ku * quadgk(@(y) atom_fun(y), -30, 30);
    exp_loss = exp_loss - atom_term;

    % Final tranche price
    price = discount * (1 - exp_loss);
end