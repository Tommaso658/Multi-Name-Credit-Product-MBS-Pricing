function [pr_HP, pr_up, pr_down, pr_HP_comp, score] = prices_vs_I_double_t( ...
                discount, Kd_calib, Ku_list, rho_LHP, recovery, p, nu_LHP, ...
                N_sim, Kd_list)
    I_vec   = round(logspace(1,3,20));
    pr_HP   = zeros(numel(I_vec),3);
    pr_up   = pr_HP;  pr_down = pr_HP;
    pr_HP_comp = tranche_prices_double_t_optimized(nu_LHP, I_vec, ...
                     Kd_calib, Kd_list, Ku_list, recovery, rho_LHP, p, discount, 0);
    for ii = 1:numel(I_vec)
        [c, cu, cd] = tCopulaCumulative(discount, Kd_calib, Ku_list, ...
                         N_sim, rho_LHP, recovery, p, I_vec(ii), nu_LHP);
        pr_HP(ii,:) = c(1:3);  pr_up(ii,:) = cu(1:3);  pr_down(ii,:) = cd(1:3);
    end
    score = score_function(pr_HP_comp(:,1:3), pr_up, pr_down);
end