function [price] = price_cumul_tranche_HP_vasicek_optimized(Kd, Ku, avg_recovery, rho, p, discount, I)
% -----------------------------------------------------------------------------------------
% price_cumul_tranche_HP_vasicek_optimized - VERSIONE OTTIMIZZATA VETTORIALE
%
% Versione ottimizzata della funzione originale con calcoli completamente vettorizzati
% per eliminare i cicli for e migliorare drasticamente le performance.
% -----------------------------------------------------------------------------------------

    % Step 1 – Normalize attachment points to loss space [0, 1]
    d = Kd / (1 - avg_recovery);
    u = Ku / (1 - avg_recovery);

    % Step 2 – Pre-compute all possible values of m (0 to I) as a vector
    m_vec = (0:I)';  % Column vector of all possible default counts
    
    % Step 3 – Pre-compute binomial coefficients for all m values
    % Uso log per evitare overflow con portfolii grandi
    log_binom_coeff = gammaln(I + 1) - gammaln(m_vec + 1) - gammaln(I - m_vec + 1);
    
    % Step 4 – Pre-compute tranche loss for all m values (vettoriale)
    loss_tranche_vec = min(max((m_vec / I - d) / (u - d), 0), 1);

    % Step 5 – Define conditional default probability p(y) 
    sqrt_rho = sqrt(rho);
    sqrt_1_minus_rho = sqrt(1 - rho);
    norminv_p = norminv(p);
    
    p_y = @(y) normcdf((norminv_p - sqrt_rho .* y) ./ sqrt_1_minus_rho);

    % Step 6 – Vectorized integrand function
    % Calcola simultaneamente per tutti i valori di m
    integrand_vectorized = @(y) sum(exp(log_binom_coeff + ...
                                      m_vec .* log(max(p_y(y), 1e-15)) + ...
                                      (I - m_vec) .* log(max(1 - p_y(y), 1e-15))) .* ...
                                      loss_tranche_vec) .* normpdf(y);

    % Step 7 – Single integration over y
    expected_loss = quadgk(integrand_vectorized, -6, 6, 'RelTol', 1e-10, 'AbsTol', 1e-12);

    % Step 8 – Final price
    price = discount * (1 - expected_loss);
end

