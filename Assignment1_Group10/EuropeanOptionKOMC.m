function  optionPrice = EuropeanOptionKOMC(F0,K,KO,B,T,sigma,N)
%Up&Out EU Call option pricing with MC approach

%INPUT PARAMETERS
% F0:    forward price
% K:     strike
% K0:    barrier
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
% N:     number of simulations

%generation of random vector from a standard gaussian distribution
g = randn(N,1);
%computation of the vector of payoffs 
XT = -sigma^2*0.5*T+sigma*sqrt(T).*g;
FT = F0 * exp(XT);
discpayoff = B * max(FT-K,0).*(FT<KO);
%computation of the price 
optionPrice = mean(discpayoff);

end

