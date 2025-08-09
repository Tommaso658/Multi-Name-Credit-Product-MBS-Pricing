function [I_vec, pr, pr_up, pr_down] = prices_vs_I_gaussian_LHP( ...
                discount, Kd_calib, Ku_list, rho_LHP, recovery, p, N_sim)
    I_vec   = round(logspace(1,3,20));
    pr      = zeros(numel(I_vec),3);
    pr_up   = pr;  pr_down = pr;
    for ii = 1:numel(I_vec)
        [c, cu, cd] = gaussian_copula(discount, Kd_calib, Ku_list, ...
                         N_sim, rho_LHP, recovery, p, I_vec(ii));
        pr(ii,:)    = c(1:3);   pr_up(ii,:) = cu(1:3);   pr_down(ii,:) = cd(1:3);
    end
end