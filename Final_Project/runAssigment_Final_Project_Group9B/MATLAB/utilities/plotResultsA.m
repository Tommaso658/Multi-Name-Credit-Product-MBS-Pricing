function plotResultsA(res)

%% Graphic defaults --------------------------------------------------------
set(groot, 'defaultFigureColor','w', ...
          'defaultAxesFontSize',12, ...
          'defaultLineLineWidth',1.5, ...
          'defaultLineMarkerSize',8);

%% A.1 – Calibration – Double t‑Student (LHP) -----------------------------
figure;
semilogy(res.pointA.nu_list, res.pointA.MSE_list_LHP,'o-'); hold on;
nu1  = res.pointA.nu_list(res.pointA.idx_min_LHP);
mse1 = res.pointA.MSE_min_LHP;
plot(nu1, mse1,'rp','MarkerSize',12,'MarkerFaceColor','r');
text(nu1+1, mse1, sprintf('$\\leftarrow\\;\\nu^* = %.2f$',nu1), ...
     'Color','r','Interpreter','latex');

xlabel('$\nu$ (degrees of freedom)','Interpreter','latex');
ylabel('MSE($\nu$)','Interpreter','latex');
title('Calibration – Double t-Student (LHP)');
grid on; grid minor;
legend({'MSE($\nu$)','$\nu^*$'}, 'Interpreter','latex', 'Location','northeast');
hold off;

%% A.2 – Improved calibration – total correlation MSE vs \nu -------------
figure;
semilogy(res.pointA.nu_imp, res.pointA.MSE_nu,'o-'); hold on;
nu2  = res.pointA.nu_star;
mse2 = res.pointA.MSE_min_imp;
plot(nu2, mse2,'rp','MarkerSize',12,'MarkerFaceColor','r');
text(nu2+1, mse2, sprintf('$\\leftarrow\\;\\nu^* = %.2f$',nu2), ...
     'Color','r','Interpreter','latex');
xlabel('$\nu$ (degrees of freedom)','Interpreter','latex');
ylabel('MSE($\rho^*,\rho^{mkt}$)','Interpreter','latex');
title('Improved calibration – total correlation MSE vs \nu', ...
    'Interpreter','tex');
grid on; grid minor;
legend({'MSE($\nu$)','$\nu^*$'}, 'Interpreter','latex', 'Location','northeast');
hold off;

%% A.3 – Comparison of the two \nu* estimation methods -------------------
figure;
semilogy(res.pointA.nu_list, res.pointA.MSE_list_LHP,'o-'); hold on;
plot(nu1, mse1,'rp','MarkerSize',12,'MarkerFaceColor','r');
plot(res.pointA.global.nu_star_global, res.pointA.global.MSE_star_global, ...
     'bs','MarkerSize',12,'MarkerFaceColor','b');

xlabel('$\nu$ (degrees of freedom)','Interpreter','latex');
ylabel('MSE','Interpreter','latex');
title('Comparison of the two \nu^* estimation methods','Interpreter', 'tex');
grid on; grid minor;
legend({'MSE curve (equity match)','$\nu^*_1$ (equity match)', ...
        '$\nu^*_2$ (global search)'}, ...
        'Interpreter','latex','Location','northeast');
hold off;

%% A.4 – Zoom on the two \nu* estimates ----------------------------------
figure;
xlo  = min(nu1, res.pointA.global.nu_star_global)-1;
xhi  = max(nu1, res.pointA.global.nu_star_global)+1;
ymin = min(mse1,res.pointA.global.MSE_star_global);

semilogy(res.pointA.nu_list, res.pointA.MSE_list_LHP,'o-'); hold on;
plot(nu1, mse1,'rp','MarkerSize',12,'MarkerFaceColor','r');
plot(res.pointA.global.nu_star_global, res.pointA.global.MSE_star_global, ...
     'bs','MarkerSize',12,'MarkerFaceColor','b');

xlim([xlo xhi]); ylim([0.8*ymin 3*ymin]);
xlabel('$\nu$ (degrees of freedom)','Interpreter','latex');
ylabel('MSE','Interpreter','latex');
title('Zoom on the two \nu^* estimates','Interpreter', 'tex');
grid on; grid minor;
legend({'MSE curve','$\nu^*_1$','$\nu^*_2$'}, ...
       'Interpreter','latex','Location','northeast');
hold off;

end
