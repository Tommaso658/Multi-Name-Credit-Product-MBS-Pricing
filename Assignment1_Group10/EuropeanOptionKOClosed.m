function optionPrice = EuropeanOptionKOClosed(F0,K,KO,B,T,sigma,flag)
%European option European Barrier up&out price with exact formula
%
%INPUT
% F0:           forward price
% r:            rate
% B:            discount factor
% K:            strike price
% KI:           Barrier strike
% T:            time-to-maturity
% sigma:        volatility

% We are seeing the Barrier Option as a bull call spread - cash or nothing
%computation of d1,d2 with strike KO
d1=log(F0/KO)/(sigma*sqrt(T))+1/2*sqrt(T)*sigma;
d2=log(F0/KO)/(sigma*sqrt(T))-1/2*sqrt(T)*sigma;
%computation of the prices 
Call_K = EuropeanOptionClosed(F0,K,B,T,sigma,1);
Call_KO = EuropeanOptionClosed(F0,KO,B,T,sigma,1);
Cash_or_nothing = B * (KO-K)*cdf('Normal',d2,0,1);

%summing the two prices
optionPrice = Call_K - Call_KO - Cash_or_nothing;

end

