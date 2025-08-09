function price = price_cumul_tranche_KL_double_t_advanced(kd, ku, recovery, rho, p, discount, I, nu, tinv_p_nu, sqrt_rho, sqrt_1_minus_rho)
% Optimized version of the pricing function with pre-computed terms

    % Normalize attachment points
    d = kd / (1 - recovery);
    u = ku / (1 - recovery);
    
    % Tranche loss profile
    loss_tr = @(z) max(0, min(1, (z - d) / (u - d)));
    
    % Conditional default probability (using pre-computed terms)
    p_y = @(y) tcdf((tinv_p_nu - sqrt_rho .* y) ./ sqrt_1_minus_rho, nu);
    
    % KL divergence
    K = @(z, y) z .* log(z ./ p_y(y)) + (1 - z) .* log((1 - z) ./ (1 - p_y(y)));
    
    % Normalization factor
    C1 = @(z) sqrt(I ./ (2 * pi .* z .* (1 - z)));
    
    % Full integrand
    integranda_zy = @(z, y) C1(z) .* exp(-I .* K(z, y)) .* loss_tr(z) .* tpdf(y, nu);
    
    % Vectorized grid computation
    dy = 0.1;
    discrete_y = (-6:dy:6)';
    y_midpoints = (discrete_y(1:end-1) + discrete_y(2:end)) / 2;
    y_intervals = discrete_y(2:end) - discrete_y(1:end-1);
    
    % Vectorized integration
    integranda_y = arrayfun(@(y_mid) quadgk(@(z) integranda_zy(z, y_mid), 0, 1), y_midpoints);
    exp_loss = dot(y_intervals, integranda_y);
    
    % Atomic correction
    atom_fun = @(y) (1 - p_y(y)).^I .* tpdf(y, nu);
    atom_term = ku * quadgk(@(y) atom_fun(y), -30, 30);
    exp_loss = exp_loss - atom_term;
    
    % Final price
    price = discount * (1 - exp_loss);
end