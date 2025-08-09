function [pr, pr_up, pr_down, pr_HP_comp, score] = prices_vs_I_gaussian_vasicek( ...
                discount, Kd_calib, Ku_list, rho_eq, recovery, p, N_sim, Kd_list)
    I_vec   = round(logspace(1,3,20));
    pr      = zeros(numel(I_vec),3);
    pr_up   = pr;  pr_down = pr;
    pr_HP_comp = tranche_prices_vasicek(I_vec, Kd_calib, Kd_list, Ku_list, ...
                       recovery, rho_eq, p, discount, 0);
    for ii = 1:numel(I_vec)
        [c, cu, cd] = gaussian_copula(discount, Kd_calib, Ku_list, ...
                         N_sim, rho_eq, recovery, p, I_vec(ii));
        pr(ii,:)      = c(1:3);
        pr_up(ii,:)   = cu(1:3);
        pr_down(ii,:) = cd(1:3);
    end
    score = score_function(pr_HP_comp(:,1:3), pr_up, pr_down);
end
