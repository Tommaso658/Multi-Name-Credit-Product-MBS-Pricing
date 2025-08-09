function gamma = vegaKO(F0,K,KO,B,T,sigma,N,flagNum)
% compute the vega of barrier option (up and in call) for CRR, MC and with the closed
% formula

%INPUT
% F0:    forward price
% K:     strike price
% KO:    barrier level
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
% N:     number of steps
% flagNum=1 CRR, flagNum=2 MC, flagNum=3 Exact

delta_sigma = 1/100; % 1%
h = 1/10; %step to compute numerically the derivative after CRR

if flagNum == 1
   %simulate prices via CRR
   Price_CRR_hplus=EuropeanOptionKOCRR(F0,K, KO,B,T,sigma+h,N);
   Price_CRR_hminus=EuropeanOptionKOCRR(F0,K, KO,B,T,sigma-h,N);
   %compute vega via centered derivatives appriximation
   gamma=((Price_CRR_hplus-Price_CRR_hminus)/(2*h))*delta_sigma;
end

h = 1/10;%step to compute numerically the derivative after MC

if flagNum == 2
        %simulate prices via MC
        price_MC_hplus = EuropeanOptionKOMC(F0,K,KO,B,T,sigma+h,N);
        price_MC_hminus = EuropeanOptionKOMC(F0,K,KO,B,T,sigma-h,N);
        %compute vega via centered derivatives approximation
        gamma = ((price_MC_hplus - price_MC_hminus)/(2*h))*delta_sigma;

end



if flagNum ==3
     %compute d1,2 for both strikes 
     d1KO=log(F0/KO)/(sigma*sqrt(T))+1/2*sqrt(T)*sigma;
     d2KO=log(F0/KO)/(sigma*sqrt(T))-1/2*sqrt(T)*sigma;
     d1K=log(F0/K)/(sigma*sqrt(T))+1/2*sqrt(T)*sigma;
     d2K=log(F0/K)/(sigma*sqrt(T))-1/2*sqrt(T)*sigma;
     %compute separately the vegas of the two call options and of the
     %digital option
     vegaCallK = B*F0*sqrt(T)*exp(-d1K^2*0.5)/sqrt(2*pi);
     vegaCallKO = B*F0*sqrt(T)*exp(-d1KO^2*0.5)/sqrt(2*pi);
     vegaDigital = B*(KO-K)*exp(-d2KO^2*0.5)/sqrt(2*pi)*(log(KO/F0)/(sigma^2*sqrt(T))-0.5*sqrt(T));
     %compute the vega of the knock-out call option seeing it as the sum of
     %the three options above 
     gamma = (vegaCallK - vegaCallKO-vegaDigital)*delta_sigma;
end



end

