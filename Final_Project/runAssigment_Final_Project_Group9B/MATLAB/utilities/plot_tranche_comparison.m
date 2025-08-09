function plot_tranche_comparison(results, I_list_E, Ku_list, tranche_prices_vasicek_LHP, tranche_prices_HP_vasicek_comparaison)
% ------------------------------------------------------------------------------------------------------------
% plot_tranche_comparison
%
% Crea un grafico di confronto per i prezzi delle tranche usando diversi modelli
% (t-Student vs Vasicek) con intervalli di confidenza e scale logaritmiche.
%
% INPUTS:
%   results                                 - struct contenente i risultati con campi:
%                                            .pointE.I_list_E, .pointE.prLi_t, 
%                                            .pointE.prLi_t_up, .pointE.prLi_t_down, etc.
%   I_list_E                               - vettore delle dimensioni del portfolio
%   Ku_list                                - vettore degli upper attachment points
%   tranche_prices_vasicek_LHP             - prezzi LHP Vasicek
%   tranche_prices_HP_vasicek_comparaison  - prezzi HP Vasicek per confronto
%
% OUTPUT:
%   Figure con 3 subplot che mostrano il confronto dei modelli
% ------------------------------------------------------------------------------------------------------------

    % Verifica che results abbia la struttura corretta
    if ~isfield(results, 'pointE')
        error('La struttura results deve contenere il campo pointE');
    end
    
    % Estrai i dati dalla struttura results
    pointE = results.pointE;
    
    % Verifica che tutti i campi necessari esistano
    required_fields = {'I_list_E', 'prLi_t', 'prLi_t_up', 'prLi_t_down', ...
                      'prTS_HP', 'prTS_up', 'prTS_down', 'prLi_vas', ...
                      'prLi_vas_up', 'prLi_vas_down', 'prTS_vas_HP'};
    
    for i = 1:length(required_fields)
        if ~isfield(pointE, required_fields{i})
            error('Campo mancante in results.pointE: %s', required_fields{i});
        end
    end
    
    % Estrai i dati
    x = I_list_E(:);
    prices_Li_HP = pointE.prLi_t;
    prices_Li_HP_down = pointE.prLi_t_down;
    prices_Li_HP_up = pointE.prLi_t_up;
    
    % Numero di tranche
    n_tranches = min(length(Ku_list), size(prices_Li_HP, 2));
    
    % Definisci i colori personalizzati per le curve
    colors_lines = [ ...
        0 0 0.6;        % Monte Carlo mean (blu scuro)
        1 0 0;          % LHP (rosso)
        0 0.6 0;        % KL (verde)
        0.3 0.75 0.93;  % KL Vasicek (blu chiaro)
        0.6 0 0.6;      % HP (magenta)
        1 0.5 0];       % HP Vasicek (arancione)
    
    % Crea la figura
    figure('Position', [100, 100, 1200, 400]);
    
    % Loop sui subplot (massimo 3 tranche)
    for k = 1:min(3, n_tranches)
        subplot(1, 3, k); 
        hold on;
        
        % Estrai i vettori per questa tranche
        y = prices_Li_HP(:, k);
        y_lower = prices_Li_HP_down(:, k);
        y_upper = prices_Li_HP_up(:, k);
        
        % Riempi l'intervallo di confidenza (zona blu)
        x_fill = [x; flipud(x)];
        y_fill = [y_upper; flipud(y_lower)];
        
        fill(x_fill, y_fill, [0.8 0.8 1], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
            'DisplayName', 'Monte Carlo CI');
        
        % Curva Monte Carlo (mean)
        plot(x, y, '-o', 'Color', colors_lines(1,:), 'LineWidth', 1.5, ...
            'MarkerSize', 4, 'DisplayName', 'Monte Carlo mean');
        
        % Linea orizzontale: LHP (fisso)
        if k <= size(tranche_prices_vasicek_LHP, 2)
            yline(tranche_prices_vasicek_LHP(1, k), '--', ...
                'Color', colors_lines(2,:), 'LineWidth', 1.2, ...
                'DisplayName', sprintf('LHP vasicek [0–%g%%]', 100 * Ku_list(k)));
        end
        
        % Curve HP
        if k <= size(tranche_prices_HP_vasicek_comparaison, 2)
            plot(x, tranche_prices_HP_vasicek_comparaison(:,k), ':', ...
                'Color', colors_lines(5,:), 'LineWidth', 1.5, ...
                'DisplayName', 'HP (t-Student)');
        end
        
        % Formatting del subplot
        title(sprintf('Tranche [0–%g%%]', 100 * Ku_list(k)), ...
            'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Number of Exposures I (log scale)', 'FontSize', 10);
        ylabel('Price (% of Notional)', 'FontSize', 10);
        
        % Scala logaritmica sull'asse x
        set(gca, 'XScale', 'log');
        
        % Griglia
        grid on;
        set(gca, 'GridAlpha', 0.3);
        
        % Legenda
        legend('Location', 'best', 'FontSize', 9);
        
        % Miglioramenti estetici
        set(gca, 'FontSize', 9);
        box on;
    end
    
    % Titolo generale
    sgtitle('Comparison of Tranche Pricing Models – t-Student vs Vasicek', ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Regola la spaziatura tra i subplot
    set(gcf, 'PaperPositionMode', 'auto');
    
    % Opzionale: salva la figura
    % saveas(gcf, 'tranche_comparison.png', 'png');
    % saveas(gcf, 'tranche_comparison.fig', 'fig');

end
