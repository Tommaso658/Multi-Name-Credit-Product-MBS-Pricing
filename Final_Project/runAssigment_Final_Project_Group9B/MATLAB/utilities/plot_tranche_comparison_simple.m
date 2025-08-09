function plot_tranche_comparison_simple(I_list, MC_mean, MC_lower, MC_upper, LHP_prices, HP_prices, Ku_list)
% ------------------------------------------------------------------------------------------------------------
% plot_tranche_comparison_simple
%
% Crea un grafico di confronto per i prezzi delle tranche usando diversi modelli
% (t-Student vs Vasicek) con intervalli di confidenza e scale logaritmiche.
%
% INPUTS:
%   I_list      - vettore delle dimensioni del portfolio (N x 1)
%   MC_mean     - prezzi Monte Carlo mean (N x n_tranches)
%   MC_lower    - intervallo di confidenza inferiore (N x n_tranches)
%   MC_upper    - intervallo di confidenza superiore (N x n_tranches)
%   LHP_prices  - prezzi LHP (valori fissi, 1 x n_tranches)
%   HP_prices   - prezzi HP t-Student (N x n_tranches)
%   Ku_list     - vettore degli upper attachment points (decimali, es. [0.03, 0.06, 0.09])
%
% OUTPUT:
%   Figure con 3 subplot che mostrano il confronto dei modelli
% ------------------------------------------------------------------------------------------------------------

    % Converti in vettore colonna se necessario
    x = I_list(:);
    
    % Numero di tranche (massimo 3 per i subplot)
    n_tranches = min(length(Ku_list), size(MC_mean, 2));
    
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
        y = MC_mean(:, k);
        y_lower = MC_lower(:, k);
        y_upper = MC_upper(:, k);
        
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
        if k <= length(LHP_prices)
            yline(LHP_prices(k), '--', ...
                'Color', colors_lines(2,:), 'LineWidth', 1.2, ...
                'DisplayName', sprintf('LHP Vasicek [0–%g%%]', 100 * Ku_list(k)));
        end
        
        % Curve HP t-Student
        if k <= size(HP_prices, 2)
            plot(x, HP_prices(:,k), ':', ...
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