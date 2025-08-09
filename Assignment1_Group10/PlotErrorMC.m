function [M, stdEstim] = PlotErrorMC(F0,K,B,T,sigma)
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
% stdEstim:  = vector of MC errors 

% Selection of M for MC method
m = 1:20;   % Row vectors
M = 2.^m;   % Row vectors with powers of 2
stdEstim = zeros(1,length(M));   % Preallocate the errors vector

for i = 1:length(M)
    
    % for each value of M, M(i) normally distributed random numbers are generated
    g = randn(M(i),1);

    % Forward price evolution at T under risk-neutral dynamics
    FT = F0*exp(-sigma^2*0.5*T+sigma*sqrt(T).*g);

    % Discounted payoff
    discpayoff = B * max(FT-K,0);
    
    % Calculates the standard error of the MC estimate
    stdEstim(i) = std(discpayoff)/sqrt(M(i));

end

% Monte Carlo Error Plot
figure()
loglog(M, stdEstim, 'o-','LineWidth',1.5,'DisplayName','MC StdError');
hold on; grid on;
xlabel('M (simulation MC)');
ylabel('Std Error');
title('Error MC vs M (log-log)');

% Let's add a reference line with slope -0.5 for show that the error
% scales ~ 1/sqrt(M)
M_mid   = M(round(end/2));   % Central point of M
err_mid = stdEstim(round(end/2));   % Error corresponding to central point of M
xRef2   = logspace(log10(min(M)), log10(max(M)), 100);    % Generate logarithmically spaced vector with 100 point
yRef2   = err_mid * (xRef2/M_mid).^(-0.5);   %straight line with slope 0.5
loglog(xRef2, yRef2, 'k--','LineWidth',1.2,'DisplayName','slope -0.5');
%loglog(M,1./sqrt(M), 'k--','LineWidth',1.2,'DisplayName','1/sqrt(M)',Color='r')   % Line sqrt(M)
legend('Location','best');

end
