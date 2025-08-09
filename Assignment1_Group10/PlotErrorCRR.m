function [M, errorCRR] = PlotErrorCRR(F0,K,B,T,sigma)
% Plot in a log-log scale the CRR error wrt M
%
% INPUT:
%   F0     = forward price at t=0 for maturity T
%   K      = strike price
%   B      = discount factor = e^{-rT}
%   T      = time to maturity
%   sigma  = volatility
%
% OUTPUT:
% M:       = vector of Time Steps 
% errorCRR = vector of CRR errors

flag = 1; %call option
% Option price with closed Formula (optionPrice = 1)
OptionPriceExact = EuropeanOptionPrice(F0,K,B,T,sigma,1,flag);

m = 1:10;   % Row vectors
M = 2.^m;   % Row vectors with powers of 2
errorCRR = zeros(1,length(m));   % Preallocate the errors vector

for i = 1:length(M)
 
    % Option price with CRR (optionPrice = 2)
    OptionPriceCRR = EuropeanOptionPrice(F0,K,B,T,sigma,2,M(i),flag);
    
    % Calculate the absolute error compared to the exact BS price
    errorCRR(i) = abs(OptionPriceCRR - OptionPriceExact);
end

% CRR Error Plot
figure()
loglog(M, errorCRR, 'o-','LineWidth',1.5,'DisplayName','CRR Error');
hold on; grid on;
xlabel('M (steps CRR)');
ylabel('|C_{CRR}(N) - C_{exact}|');
title('Error CRR vs M (log-log)');

% Let's add a reference line with slope -1 for show that the error
% scales ~ 1/M
M_mid   = M(round(end/2));   % Central point of M
err_mid = errorCRR(round(end/2));   % Error corresponding to central point of M
xRef    = logspace(log10(min(M)), log10(max(M)), 100);   % Generate logarithmically spaced vector with 100 point
yRef    = err_mid * (xRef/M_mid).^(-1);   % Straight line with slope 1
loglog(xRef, yRef, 'k--','LineWidth',1.2,'DisplayName','slope -1');
%loglog(M,1./M, 'k--','LineWidth',1.2,'DisplayName','1/M',Color='r');   % Line 1/sqrt(M)
legend('Location','best');

end

