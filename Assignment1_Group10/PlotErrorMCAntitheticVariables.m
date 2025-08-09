function [M, stdEstimAV] = PlotErrorMCAntitheticVariables(F0,K,B,T,sigma)
% Plot in a log-log scale the MC error wrt M
%
% INPUT:
%   F0     = forward price at t=0 for maturity T
%   K      = strike price
%   B      = discount factor = e^{-rT}
%   T      = time to maturity
%   sigma  = volatility
%
% OUTPUT:
% M          = vector of Time Steps 
% stdEstimAV = vector of MC errors using AV

% Selection of M for MC method
m = 1:20;   % Row vectors
M = 2.^m;   % Row vectors with powers of 2
% stdEstim = zeros(1,length(M));   % Preallocate the errors vector
stdEstimAV = zeros(1,length(m));   % Preallocate the AV errors vector

for i = 1:length(m)
    
    % for each value of M, M(i) normally distributed random numbers are generated
    g1 = randn(M(i),1);
    g2 = -g1;   %Antithetic variables as oppoisite of g

    % Forward price evolution at T under risk-neutral dynamics
    FT1 = F0*exp(-sigma^2*0.5*T+sigma*sqrt(T).*g1);
    FT2 = F0*exp(-sigma^2*0.5*T+sigma*sqrt(T).*g2);

    % Discounted payoff
    discpayoff1 = B * max(FT1-K,0);
    discpayoff2 = B * max(FT2-K,0);
    % Payoff of AV
    discpayoff = (discpayoff1+discpayoff2)/2;

    % Compute the standard error of the MC estimate
    stdEstimAV(i) = std(discpayoff)/sqrt(M(i));
    %stdEstim(i) = std(discpayoff1)/sqrt(M(i));

end

end
