function tranche_prices = tranche_prices_double_t_optimized(nu_optimal, I_list, Kd_calibration, Kd_list, Ku_list, recovery, corr_model, p, discount, flag)

% -------------------------------------------------------------------------------------------------------
% tranche_prices_double_t_optimized - VERSIONE OTTIMIZZATA VETTORIALE
%
% Versione ottimizzata che elimina i cicli for annidati usando calcoli vettoriali.
% -------------------------------------------------------------------------------------------------------

    % Pre-allocate result matrix
    tranche_prices = zeros(length(I_list), length(Ku_list));
    
    % Round I_list once
    I_list_rounded = round(I_list);

    if flag == 0
        % Exact pricing - Vectorized over tranches for each I
        for j = 1:length(I_list)
            I_val = I_list_rounded(j);
            
            % Vectorized computation over all Ku values
            tranche_prices(j, :) = arrayfun(@(ku) ...
                price_cumul_tranche_HP_double_t_student_optimized( ...
                    Kd_calibration, ku, recovery, corr_model, p, discount, I_val, nu_optimal), ...
                Ku_list);
        end
        
    elseif flag == 1
        % KL Approximation - Vectorized over tranches for each I
        for j = 1:length(I_list)
            I_val = I_list_rounded(j);
            
            % Vectorized computation over all Ku values
            tranche_prices(j, :) = arrayfun(@(ku) ...
                price_cumul_tranche_KL_double_t_optimized( ...
                    Kd_calibration, ku, recovery, corr_model, p, discount, I_val, nu_optimal), ...
                Ku_list);
        end
        
    elseif flag == 2
        % LHP Approximation - VERSIONE COMPLETAMENTE OTTIMIZZATA
        % Calcola tutti i tranches simultaneamente con una singola chiamata
        
        % Calcola i prezzi per tutti i tranches in una volta sola
        lhp_prices = price_LHP_t_student_optimized(discount, recovery, Ku_list, ...
                                                  Kd_calibration, corr_model, p, nu_optimal);
        
        % Replica i risultati per tutti i valori di I (poiché LHP non dipende da I)
        tranche_prices = repmat(lhp_prices, length(I_list), 1);
    end
end