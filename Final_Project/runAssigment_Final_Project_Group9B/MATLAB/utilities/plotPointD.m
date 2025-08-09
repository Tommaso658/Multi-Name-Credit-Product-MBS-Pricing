function plotPointD(res)
% Point D (Vasicek) + Point E (Li & t‑Student comparisons)
setGraphicDefaults;
set(groot,'defaultTextInterpreter','latex');   % (una sola volta, all’inizio)

%% D.1 – Vasicek bar chart (I = 500)
figure;
barData = [res.pointD.prVas_HP(1:3); res.pointD.prVas_KL(1:3); res.pointD.prVas_LHP(1:3)];
bar(barData','grouped'); hold on;
set(gca,'XTickLabel',{'0–3%','0–6%','0–9%'});
leg = {'Exact','KL','LHP'}; cols=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125];
for i=1:3, b=findobj(gca,'Type','Bar','SeriesIndex',i); set(b,'FaceColor',cols(4-i,:)); end %#ok<*AGROW>
xlabel('Tranche'); ylabel('Price (% Notional)'); grid on; grid minor;
title(sprintf('Vasicek Model – \\rho_{equity}=%.3f, I=%d', ...
    res.pointD.prVas_HP(1), res.pointD.I_vas), ...
    'Interpreter','tex');
legend(fliplr(leg),'Location','northwest');

%% D.2 – Equity tranche t‑Student vs Vasicek
figure; hold on;
semilogx(res.pointD.I_list, res.pointB.pricesHP(:,1),'o-','Color',[0 0.447 0.741],'DisplayName','t‑Student (Exact)');
semilogx(res.pointD.I_list, res.pointD.prVas_HP_list(:,1),'s-','Color',[0.85 0.325 0.098],'DisplayName','Vasicek (KL)');
yline(res.pointB.pricesLHP(1,1),'--','Color','k','LineWidth',1.5,'DisplayName','LHP');
ax=gca; ax.XScale='linear'; ax.XTick=[200 400 600 800 1000]; ax.XLim=[0 1000]; ax.YTick=0.30:0.05:0.60; ax.YLim=[0.30 0.60];
xlabel('Number of Exposures (I)'); ylabel('Equity Tranche Price (% Notional)');
title('Equity Tranche: t‑Student vs Vasicek','Interpreter','latex');
legend('Location','northeast'); grid on; grid minor;