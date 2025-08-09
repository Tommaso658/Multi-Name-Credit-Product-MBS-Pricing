
function price = price_cumul_tranche_KL_vasicek_optimized(kd, ku, recovery, rho, p, disc, I)
% -----------------------------------------------------------------------------------------
% price_cumul_tranche_KL_vasicek_optimized - VERSIONE OTTIMIZZATA VETTORIALE
%
% Versione ottimizzata che usa calcoli vettoriali per l'integrazione su griglia.
% -----------------------------------------------------------------------------------------

    % Transform attachment points
    d = kd / (1 - recovery);
    u = ku / (1 - recovery);

    % Loss profile function
    loss_tr = @(z) max(0, min(1, (z - d) / (u - d)));

    % Pre-compute constants
    sqrt_rho = sqrt(rho);
    sqrt_1_minus_rho = sqrt(1 - rho);
    norminv_p = norminv(p);
    
    % Conditional default probability
    p_y = @(y) normcdf((norminv_p - sqrt_rho .* y) ./ sqrt_1_minus_rho);

    % KL divergence and normalization
    K = @(z, y) z .* log(z ./ p_y(y)) + (1 - z) .* log((1 - z) ./ (1 - p_y(y)));
    C1 = @(z) sqrt(I ./ (2 * pi .* z .* (1 - z)));

    % Vectorized integrand
    integranda_zy = @(z, y) C1(z) .* exp(-I .* K(z, y)) .* loss_tr(z) .* normpdf(y);

    % Optimized grid integration
    dy = 0.1;
    y_grid = (-6:dy:6)';
    y_mid = (y_grid(1:end-1) + y_grid(2:end)) / 2;  % Midpoints vectorized
    dy_vec = diff(y_grid);  % Differences vectorized
    
    % Vectorized integration over z for all y midpoints
    integranda_y = arrayfun(@(y_val) quadgk(@(z) integranda_zy(z, y_val), 0, 1, ...
                                          'RelTol', 1e-8), y_mid);
    
    % Vectorized sum for expected loss
    exp_loss = sum(dy_vec .* integranda_y);

    % Atomic correction
    atom_fun = @(y) (1 - p_y(y)).^I .* normpdf(y);
    atom_term = ku * integral(atom_fun, -Inf, Inf, 'RelTol', 1e-8);
    exp_loss = exp_loss - atom_term;

    % Final price
    price = disc * (1 - exp_loss);
end
