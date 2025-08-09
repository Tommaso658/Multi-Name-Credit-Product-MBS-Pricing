function OptionPrice = EuropeanOptionCRR(F0, K, B, T, sigma, N, flag)

% EuropeanOptionCRR - Price of a European call/put option (flag=1/-1)
% using the CRR (Cox-Ross-Rubinstein) binomial model with N steps.
%
% INPUT:
%   F0     = forward price at t=0 for maturity T
%   K      = strike price
%   B      = discount factor = e^{-rT}
%   T      = time to maturity
%   sigma  = volatility
%   N      = number of steps in the binomial tree
%   flag   = 1 for call, -1 for put
%
% OUTPUT:
%   OptionPrice = value of the option at t=0


    % Time length of each step
    deltaT = T / N;
    
    % Up and down factors
    u = exp(sigma * sqrt(deltaT));
    d = 1 / u;
    
    % Risk-neutral probability consistent with a forward that has zero drift
    p = (1 - d) / (u - d);

    % Payoffs at the final nodes
    payoffs = zeros(N+1, 1);
    for k = 0:N
        % Number of up-moves = k, down-moves = N-k
        % Final forward price at this node
        F_final = F0 * (u^k) * (d^(N-k));
        
        % Compute payoff based on call or put
        if flag == 1
            % Call
            payoffs(k+1) = max(F_final - K, 0);
        else
            % Put
            payoffs(k+1) = max(K - F_final, 0);
        end
    end

    % Binomial risk-neutral expected value of the payoff:
    % sum_{k=0..N} [C(N,k) * p^k * (1-p)^(N-k) * payoff(k)]
    OptionPrice = 0;
    for k = 0:N
        OptionPrice = OptionPrice + nchoosek(N, k) * p^k * (1 - p)^(N - k) * payoffs(k+1);
    end

    % % Discount from time T to 0 using B = e^{-rT}
    OptionPrice = B * OptionPrice;
end
