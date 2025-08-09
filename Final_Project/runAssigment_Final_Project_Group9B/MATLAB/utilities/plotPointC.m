function plotPointC(res)
% Point C – KL‑based calibration MSE vs nu
setGraphicDefaults;

figure;
semilogy(res.pointC.nu_list, res.pointC.MSE_list_KL,'o-'); hold on;
nu3  = res.pointC.nu_list(res.pointC.idx_min_KL);
mse3 = res.pointC.MSE_min_KL;
plot(nu3,mse3,'rp','MarkerSize',12,'MarkerFaceColor','r');
text(nu3+1,mse3,sprintf('  $\\nu^* = %.2f$',nu3), 'Color','r','Interpreter','latex');
xlabel('$\nu$ (degrees of freedom)','Interpreter','latex');
ylabel('MSE','Interpreter','latex');
title('Calibration – Double t‑Student (KL approximation)','Interpreter','latex');
legend({'MSE($\nu$)','$\nu^*$'},'Interpreter','latex','Location','northeast');
grid on; grid minor;
end