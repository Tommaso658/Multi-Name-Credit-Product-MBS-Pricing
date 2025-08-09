function [price] = price_cumul_tranche_HP_double_t_student_optimized(Kd, Ku, avg_recovery, rho, p, discount, I, nu)
% -----------------------------------------------------------------------------------------
% price_cumul_tranche_HP_double_t_student_optimized - VERSIONE OTTIMIZZATA VETTORIALE
%
% Versione ottimizzata della funzione originale con calcoli completamente vettorizzati.
% -----------------------------------------------------------------------------------------

    % Step 1 – Normalize the tranche bounds
    d = Kd / (1 - avg_recovery);
    u = Ku / (1 - avg_recovery);

    % Step 2 – Pre-compute vectors for all m values
    m_vec = (0:I)';
    
    % Step 3 – Pre-compute log binomial coefficients
    log_binom_coeff = gammaln(I + 1) - gammaln(m_vec + 1) - gammaln(I - m_vec + 1);
    
    % Step 4 – Pre-compute tranche losses (vettoriale)
    loss_tranche_vec = min(max((m_vec / I - d) / (u - d), 0), 1);

    % Step 5 – Pre-compute constants for t-distribution
    sqrt_rho = sqrt(rho);
    sqrt_1_minus_rho = sqrt(1 - rho);
    tinv_p = tinv(p, nu);
    
    % Conditional default probability
    p_y = @(y) tcdf((tinv_p - sqrt_rho .* y) ./ sqrt_1_minus_rho, nu);

    % Step 6 – Vectorized integrand
    integrand_vectorized = @(y) sum(exp(log_binom_coeff + ...
                                      m_vec .* log(max(p_y(y), 1e-15)) + ...
                                      (I - m_vec) .* log(max(1 - p_y(y), 1e-15))) .* ...
                                      loss_tranche_vec) .* tpdf(y, nu);

    % Step 7 – Integration
    expected_loss = quadgk(integrand_vectorized, -6, 6, 'RelTol', 1e-10, 'AbsTol', 1e-12);

    % Step 8 – Final price
    price = discount * (1 - expected_loss);
end
