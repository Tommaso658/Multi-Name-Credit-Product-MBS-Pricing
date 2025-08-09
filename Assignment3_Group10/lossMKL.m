function integr = lossMKL(I, pi, Ku, Kd, p, rho)
    % Compute Loss function by kullback leiber integration
    % I   : total obligors in the reference portfolio
    % pi  : recovery of obligor's debt
    % Ku  : amount of losses a portfolio can suffer before tranche's is
    %       wiped out (TRANCHE SENZA VALORE)
    % Kd  : amount of losses a portfolio can suffer before tranche's
    %       notional is eroded (TRANCHE INIZIA A PERDERE)
    % p   : defoult probability of each mortgage
    % rho : correlation betweem mortage


    % We define the normalization functions
    C1 = @(z) sqrt(I*(2 * pi * (1 - z) .* z).^-1); 
    cond_prob = @(y) normcdf((norminv(p) - sqrt(rho) * y) ./ sqrt(1 - rho)); 
    K = @(z, y) z.*log(z/cond_prob(y)) + (1 - z).*log((1 - z)/(1 - cond_prob(y))); 
    
    
    % We use Cavalieri Simpson's method to define D
    H = 0.005;
    N = 1 / H;
    z = linspace(0.00001, 0.99999, N);
    
    % Define D as a function of y using cavalieri simpson
    z_mid = (z(1:end-1)+z(2:end))/2;
    D = @(y) sum(C1(z(1:end-1)).*exp(-I*K(z(1:end-1),y))+C1(z(2:end)).*exp(-I*K(z(2:end),y))+4*C1(z_mid).*exp(-I*K(z_mid,y)))*H/6;
    
    % Define C(y, z)
    C = @(z, y) C1(z)/D(y); 
    

    % Define the loss as a function of z
    u = Ku/(1-pi);
    d = Kd/(1-pi);
    L = @(z) (min(max(z - d, 0), u - d)) / (u - d); 
    
    
    % Define the integral in the expected value function of LCK as a function of y
    LCK = @(z, y) L(z) .* C(z, y) .* exp(-I * K(z, y)); 
    
    
    % Define the integral of LCK using simpons
    integralLCK = @(y) sum(LCK(z(1:end-1), y) + LCK(z(2:end), y) + ...
                        4 * LCK(z(1:end-1)/2 + z(2:end)/ 2, y)) * H / 6;
    
    %we now compute the final integral using again cavalieri simpson
    H = 0.005;
    N = 12/H;
    y = linspace(-6, 6, N);
    cavSim = 0;
    for index = 1:N-1
        addend1 = integralLCK(y(index))*normpdf(y(index));
        addend2 = integralLCK(y(index+1))*normpdf(y(index));
        addend3 = 4*integralLCK(y(index+1)/2+y(index)/2)*normpdf(y(index+1)/2+y(index)/2);
        cavSim = cavSim +addend1+addend2+addend3;
    end 
    integr = cavSim*H/6;

end