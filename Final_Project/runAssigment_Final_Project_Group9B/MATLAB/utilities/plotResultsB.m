function plotResultsB(res)
% Point B figures (B.2 & B.5) – still in percent of notional
setGraphicDefaults;

upp = [3 6 9];
I   = res.pointB.I_list;
HP  = res.pointB.pricesHP;
KL  = res.pointB.pricesKL;
LHP = res.pointB.pricesLHP(1,:);
nu_star = res.pointA.nu_list(res.pointA.idx_min_LHP);

%% B.2 – Convergence, one subplot per tranche
figure;
for k = 1:3
    subplot(1,3,k); hold on;
    semilogx(I, HP(:,k),'o-','Color',[0 0.447 0.741],'DisplayName','Exact');
    semilogx(I, KL(:,k),'s-.','Color',[0.85 0.325 0.098],'DisplayName','KL');
    yline(LHP(k),'--','Color','k','LineWidth',1.5,'DisplayName','LHP');
    title(sprintf('Tranche 0–%d%%',upp(k)));
    xlabel('Number of Exposures (I)');
    if k==1, ylabel('Tranche Price (% Notional)'); end
    grid on; grid minor;
    ax = gca; ax.XTick=[10 50 100 200 400 600 800 1000]; ax.YTick=0.30:0.05:0.65;
    legend('Location','northeast');
end
sgtitle(sprintf('Convergence to LHP (t‑copula, \\nu = %.2f)',nu_star),'Interpreter','tex');

%% B.5 – Overlay, linear X‑scale
figure; hold on;
styles = {'o-','s-','^-'}; cols=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125];
for k = 1:3
    plot(I, HP(:,k), styles{k}, 'Color',cols(k,:), 'MarkerFaceColor','none', ...
         'DisplayName',sprintf('Exact 0–%d%%',upp(k)));
    plot(I, KL(:,k), styles{k}, 'LineStyle','-.','Color',cols(k,:), ...
         'MarkerFaceColor','none','DisplayName',sprintf('KL    0–%d%%',upp(k)));
    yline(LHP(k),'--','Color',cols(k,:),'LineWidth',1.5,'DisplayName',sprintf('LHP   0–%d%%',upp(k)));
end
ax = gca; ax.XScale='linear'; ax.XLim=[min(I) max(I)]; ax.XTick=[200 400 600 800 1000];
ax.YLim=[0.30 0.65]; ax.YTick=0.30:0.05:0.65;
xlabel('Number of Exposures (I)'); ylabel('Tranche Price (% Notional)');
title('Convergence of Exact and KL Prices to LHP (t‑copula, \nu^*)','Interpreter','tex');
legend('Location','northeast'); grid on; grid minor;
end