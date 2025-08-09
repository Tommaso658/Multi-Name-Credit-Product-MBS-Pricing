function [MSE_list, rho_list, idx_min] = calibrate_phase(nu_list, discounts, corr_mkt, Kd_calibration, Ku_list, recovery, p, I, displayFlag, is_coarse_phase, initial_rho_guess)
    % Funzione ausiliaria per calibrare una fase del grid search
    
    if nargin < 11
        initial_rho_guess = 0.3; % Default guess
    end
    
    MSE_list = inf(1, length(nu_list));
    rho_list = zeros(1, length(nu_list));
    
    % Pre-calcolo del prezzo equity di mercato (fatto una sola volta)
    price_equity_mkt = price_cumul_tranche_KL_vasicek(Kd_calibration, Ku_list(1), ...
                                                      recovery, corr_mkt(1), ...
                                                      p, discounts, I);
    
    % Pre-calcolo dei prezzi di mercato per tutte le tranche (fatto una sola volta)
    market_prices = zeros(1, length(Ku_list));
    for j = 1:length(Ku_list)
        market_prices(j) = price_cumul_tranche_KL_vasicek(Kd_calibration, Ku_list(j), ...
                                                          recovery, corr_mkt(j), p, ...
                                                          discounts, I);
    end
    
    % Setup solver options con tolleranze adattive
    if is_coarse_phase
        % Tolleranze più rilassate per la fase grossolana
        tol_fun = 1e-8;
        max_iter = 200;
    else
        % Tolleranze più strette per la fase fine
        tol_fun = 1e-12;
        max_iter = 500;
    end
    
    options = optimoptions('fsolve', ...
                           'TolFun', tol_fun, ...
                           'Display', displayFlag, ...
                           'MaxIterations', max_iter, ...
                           'Algorithm', 'levenberg-marquardt'); % Spesso più veloce
    
    % Variabili per early stopping
    consecutive_no_improvement = 0;
    best_mse_so_far = inf;
    EARLY_STOP_THRESHOLD = 5; % Ferma dopo 5 iterazioni senza miglioramento
    
    current_rho_guess = initial_rho_guess;
    
    for i = 1:length(nu_list)
        nu = nu_list(i);
        
        % Early stopping per la fase grossolana
        if is_coarse_phase && consecutive_no_improvement >= EARLY_STOP_THRESHOLD
            break;
        end
        
        % Definisci l'equazione da risolvere
        equation = @(rho) price_cumul_tranche_KL_double_t(Kd_calibration, Ku_list(1), ...
                                                          recovery, rho, p, ...
                                                          discounts, I, nu) ...
                          - price_equity_mkt;
        
        % Risolvi per rho
        try
            rho_star = fsolve(equation, current_rho_guess, options);
            
            % Aggiorna il guess per la prossima iterazione (warm start)
            current_rho_guess = rho_star;
            
        catch
            if ~strcmp(displayFlag, 'off')
                warning('Convergence problem for nu = %.2f', nu);
            end
            consecutive_no_improvement = consecutive_no_improvement + 1;
            continue;
        end
        
        % Calcola MSE in modo ottimizzato
        mse = compute_mse_optimized(Kd_calibration, Ku_list, recovery, rho_star, p, discounts, I, nu, market_prices);
        
        MSE_list(i) = mse;
        rho_list(i) = rho_star;
        
        % Controlla per early stopping
        if mse < best_mse_so_far
            best_mse_so_far = mse;
            consecutive_no_improvement = 0;
        else
            consecutive_no_improvement = consecutive_no_improvement + 1;
        end
    end
    
    [~, idx_min] = min(MSE_list);
end
