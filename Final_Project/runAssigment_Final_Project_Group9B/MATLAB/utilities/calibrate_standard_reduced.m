function [MSE_min, rho_model, idx_min, nu_list, MSE_list] = calibrate_standard_reduced(nu_list, discounts, corr_mkt, Kd_calibration, Ku_list, recovery, p, I, displayFlag)
    % Fallback al metodo standard ma con grid ridotto
    
    MSE_list = inf(1, length(nu_list));
    MSE_min = inf;
    idx_min = 1;
    rho_model = 0.3;
    
    price_equity_mkt = price_cumul_tranche_KL_vasicek(Kd_calibration, Ku_list(1), ...
                                                      recovery, corr_mkt(1), ...
                                                      p, discounts, I);
    
    options = optimoptions('fsolve', ...
                           'TolFun', 1e-10, ...
                           'Display', displayFlag, ...
                           'MaxIterations', 400);
    
    for i = 1:length(nu_list)
        nu = nu_list(i);
        
        equation = @(rho) price_cumul_tranche_KL_double_t(Kd_calibration, Ku_list(1), ...
                                                          recovery, rho, p, ...
                                                          discounts, I, nu) ...
                          - price_equity_mkt;
        
        try
            rho_star = fsolve(equation, 0.3, options);
        catch
            continue;
        end
        
        mse = 0;
        for j = 2:length(Ku_list)
            price_t = price_cumul_tranche_KL_double_t(Kd_calibration, Ku_list(j), ...
                                                      recovery, rho_star, p, ...
                                                      discounts, I, nu);
            price_v = price_cumul_tranche_KL_vasicek(Kd_calibration, Ku_list(j), ...
                                                     recovery, corr_mkt(j), p, ...
                                                     discounts, I);
            mse = mse + (price_t - price_v)^2;
        end
        
        MSE_list(i) = mse;
        
        if mse < MSE_min
            MSE_min = mse;
            idx_min = i;
            rho_model = rho_star;
        end
    end
end