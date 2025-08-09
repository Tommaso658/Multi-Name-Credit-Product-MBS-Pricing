function OptionPrice = EuropeanOptionMC(F0, K, B, T, sigma, M, flag)

% EuropeanOptionMC - Price of a European call/put option (flag=1/-1)
% using a Monte Carlo simulation with M runs.
%
% INPUT:
%   F0     = forward price at t=0 for maturity T
%   K      = strike price
%   B      = discount factor = e^{-rT}
%   T      = time to maturity
%   sigma  = volatility
%   M      = number of Monte Carlo simulations
%   flag   = 1 for call, -1 for put
%
% OUTPUT:
%   OptionPrice = value of the option at t=0

    % Preallocate the payoff vector
    payoffs = zeros(M, 1);

    for i = 1:M
        % Random draw Z ~ N(0,1)
        Z = randn;

        % Price evolution at T under risk-neutral dynamics
        S_T = F0 * exp(-0.5 * sigma^2 * T + sigma * sqrt(T) * Z);

        % Compute payoff
        if flag == 1
            % Call
            payoffs(i) = B*max(S_T - K, 0);
        else
            % Put
            payoffs(i) = B*max(K - S_T, 0);
        end
    end

    % Avarege of payoffs
    OptionPrice = mean(payoffs);
end
