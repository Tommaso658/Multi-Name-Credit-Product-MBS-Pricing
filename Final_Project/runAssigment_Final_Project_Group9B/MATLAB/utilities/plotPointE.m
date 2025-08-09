function plotPointE(res)

setGraphicDefaults;

%% E – Plot: Li's Gaussian Copula – Prices vs I (t-Student Corr) 1
figure; hold on;
styles_E = {'-o','-s','-^'};            % marker per le 3 tranche
Ku_pct   = [3,6,9];                    % percentuali superior attachment
LHP      = res.pointB.pricesLHP(1,:);  % LHP preso da Point B

for k = 1:3
    % Copia prezzi dal struct
    semilogx(res.pointE.I_list_E, res.pointE.prLi_t(:,k), styles_E{k}, ...
             'LineWidth', 1.1, ...
             'DisplayName', sprintf('Li–G 0–%d%%', Ku_pct(k)));
    % Linea orizzontale LHP
    yline(LHP(k), '--', 'LineWidth', 1.1, ...
          'DisplayName', sprintf('LHP 0–%d%%',      Ku_pct(k)));
end

grid on; grid minor;
xlabel('Number of Exposures I (log scale)');
ylabel('Tranche Price (% of Notional)');
title('Li Model (Gaussian Copula) – Prices vs I (t-Student Corr)');
legend('Location','northeast');
hold off;

%% E – Plot: t_student Copula – Prices vs I - t-student Correlation 2
I       = res.pointE.I_list_E;         % vettore delle exposure
P_MC    = res.pointE.prTS_HP;          % Monte Carlo mean (t-Student)
P_HP    = res.pointE.pr_comp_ts;      % HP deterministico t-Student
P_LHP   = res.pointB.pricesLHP;        % LHP orizzontale
Ku_pct  = [3,6,9];                     % attachement points in %

figure; hold on;

% colori originali (lines colormap)
cols = lines(3);
% marker originali
mk  = {'o','s','^'};

for k = 1:3
    % 1) Monte Carlo mean
    semilogx( I, P_MC(:,k), ...
        'LineStyle','-', ...
        'Marker',    mk{k}, ...
        'MarkerFaceColor','w', ...
        'Color',     cols(k,:), ...
        'LineWidth', 1.2, ...
        'DisplayName', sprintf('Copula [0–%d%%]',Ku_pct(k)) );

    % 2) HP deterministico t-Student
    semilogx( I, P_HP(:,k), ...
        'LineStyle','--', ...
        'Marker',    mk{k}, ...
        'Color',     cols(k,:), ...
        'LineWidth', 1.2, ...
        'DisplayName', sprintf('HP [0–%d%%]',Ku_pct(k)) );

    % 3) LHP orizzontale
    yline( P_LHP(1,k), ...
        'LineStyle','--', ...
        'Color',     cols(k,:), ...
        'LineWidth', 1.2, ...
        'DisplayName', sprintf('LHP [0–%d%%]',Ku_pct(k)) );
end

grid on;
xlabel('Number of Exposures I (log scale)');
ylabel('Tranche Price (% of Notional)');
title('Double t-Student Model – Prices comparison vs I');

% xtick e ytick fissi come nello screenshot
set(gca, 'XTick', [10 50 100 200 400 600 800 1000]);
set(gca, 'YTick', 0.30:0.05:0.65);

xlim([min(I) max(I)]);
legend('Location','northeast');
hold off;

%% E – Plot: Gaussian Copula – Prices vs I - vasicek equity Correlation 3
I      = res.pointE.I_list_E;           % rimane invariato
P_HP   = res.pointE.prLi_vas;            % Monte Carlo mean
P_HPcmp= res.pointE.prTS_vas_HP;        % HP vasicek (deterministico)

P_LHP = res.pointD.prVas_LHP(1, :)';   % <-- campo giusto

Ku = [0.03,0.06,0.09]*100;                   % attachment in % per i tre tranches

figure; hold on;
styles = {'-o','-s','-^'};                   % equity, mezz1, mezz2
colors = lines(3);                          % stesse palette di prima

for k = 1:3
    % Monte Carlo mean (Gaussian)
    semilogx(I, P_HP(:,k), styles{k}, ...
             'Color', colors(k,:), ...
             'LineWidth',1.2, ...
             'MarkerFaceColor','w', ...
             'DisplayName', sprintf('Monte Carlo mean [0–%d%%]',Ku(k)));

    % HP double t-Student deterministico
    % Estrai marker e definisci stile tratteggiato
    marker = styles{k}(2);       % per '-o' => 'o', per '-s' => 's', ecc.
    lineStyle = '--';
    
    semilogx(I, P_HPcmp(:,k), ...
             'LineStyle', lineStyle, ...
             'Marker',   marker, ...
             'Color',    colors(k,:), ...
             'LineWidth',1.2, ...
             'DisplayName', sprintf('HP (Vasicek) [0–%d%%]',Ku(k)));

    % LHP tratteggiata
    yline(P_LHP(k), '--', ...               % <-- usa P_LHP(k), non P_LHP(1,k)
          'Color', colors(k,:), ...
          'LineWidth',1.2, ...
          'DisplayName', sprintf('LHP [0–%d%%]',Ku(k)));
end

grid on;
xlabel('Number of Exposures I (log scale)');
ylabel('Tranche Price (% of Notional)');
title('Li Model (Gaussian Copula) – Prices vs I');
legend('Location','northeast');
xlim([min(I) max(I)]);
set(gca,'XTick',[10 50 100 200 400 600 800 1000]);
set(gca,'YTick',0.30:0.05:0.65);
hold off;

%% ================= Point E – Comparison Models ================
% Fig.13 – t-Student Copula – Prices vs I (con CI shaded)
figure;
cols = lines(3);
alpha_fill = 0.85;                         % più marcata

for k = 1:3
    subplot(1,3,k); hold on;
    % Estrazione serie
    xv    = I(:);
    yup   = res.pointE.prTS_up(:,k);
    ylow  = res.pointE.prTS_down(:,k);
    ymean = res.pointE.prTS_HP(:,k);
    lhp   = LHP(k);
    yhp = res.pointE.pr_comp_ts(:,k);

    % area di confidenza (blu scuro e più opaco)
    fill([xv; flipud(xv)], [yup; flipud(ylow)], ...
         [0.4 0.6 1], 'EdgeColor','none','FaceAlpha', alpha_fill);

    % media MC t-Student
    plot(xv, ymean, '-o', ...
         'Color', cols(k,:), 'LineWidth',1.5, ...
         'MarkerFaceColor','w');

    % HP t_student
    plot(xv, yhp, '-s', ...
         'Color', [0.7 0 0], 'LineWidth', 1.8, ...
         'MarkerFaceColor','w');

    % linea LHP
    yline(lhp, '--', ...
          'Color', cols(k,:), 'LineWidth',1.5);

    % formattazione
    title(sprintf('Tranche [0–%d%%]', Ku_pct(k)));
    xlabel('I (log scale)'); ylabel('Price (% Notional)');
    set(gca,'XScale','log', ...
             'XTick',[10 50 100 200 400 600 800 1000], ...
             'YTick',0.30:0.05:0.65, ...
             'XLim',[min(xv) max(xv)]);
    grid on; grid minor;

    % legenda solo nel primo subplot
    %if k==1
      legend({'CI 95%', 'Monte Carlo mean', 'HP t-Student', ...
                sprintf('LHP [0–%d%%]', Ku_pct(k))}, ...
                'Location', 'northeast');
    %end
end
sgtitle('Comparison of Tranche Pricing Models – t-Student Copula vs HP Double t-Student');


end