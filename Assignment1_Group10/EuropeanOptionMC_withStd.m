function [OptionPrice, stdError] = EuropeanOptionMC_withStd(F0, K, B, T, sigma, M, flag)
% EuropeanOptionMC_withStd - Prezzo di un'opzione europea call/put 
% (flag=1/-1) con simulazione Monte Carlo a M simulazioni.
% Restituisce anche lo standard error della stima.

    if nargin < 7
        flag = 1; % default: call
    end

    payoffs = zeros(M,1);

    for i = 1:M
        Z = randn;  
        % Evoluzione forward-based
        S_T = F0 * exp( -0.5*sigma^2*T + sigma*sqrt(T)*Z );

        if flag == 1
            payoffs(i) = max(S_T - K, 0);
        else
            payoffs(i) = max(K - S_T, 0);
        end
    end

    % Prezzo = sconto * media dei payoff
    OptionPrice = B * mean(payoffs);

    % Standard error = B * (std payoff / sqrt(M))
    payoffStd = std(payoffs); 

    stdError  = B * payoffStd / sqrt(M);
end
