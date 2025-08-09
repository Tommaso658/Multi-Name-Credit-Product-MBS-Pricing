function prices = compute_all_tranche_prices(kd, ku_list, recovery, rho, p, discounts, I, nu)
% Vectorized computation of prices for multiple tranches
% This function can be further optimized by pre-computing common terms

    n_tranches = length(ku_list);
    prices = zeros(1, n_tranches);
    
    % Pre-compute common terms that don't depend on ku
    tinv_p_nu = tinv(p, nu);
    sqrt_rho = sqrt(rho);
    sqrt_1_minus_rho = sqrt(1 - rho);
    
    % For each tranche, compute price using vectorized operations where possible
    for j = 1:n_tranches
        prices(j) = price_cumul_tranche_KL_double_t_advanced(kd, ku_list(j), recovery, ...
                                                             rho, p, discounts, I, nu, ...
                                                             tinv_p_nu, sqrt_rho, sqrt_1_minus_rho);
    end
end

