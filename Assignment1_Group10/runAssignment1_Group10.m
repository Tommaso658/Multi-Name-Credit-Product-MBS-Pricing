% Assignment_1
%  Group 10, AA2024-2025
%
%  TBM (To Be Modified): Modify & Add where needed
clear all
close all 
clc
warning off

%% Pricing parameters
S0 = 1;             % Initial underlying price (S0 = 1 Euro)
K = 1.05;           % Strike price
r = 0.025;          % Risk-free interest rate (annualized)
TTM = 1/3;          % Time to maturity (in years) = 4 months
sigma = 0.21;       % Volatility (annualized)
flag = 1;           % Option type flag: 1 for Call, -1 for Put
d = 0.02;           

% Set the seed for Monte Carlo simulations to ensure reproducibility
rng(10);

%% Quantity of interest
B=exp(-r*TTM); % Discount

%% a) Pricing 
F0=S0*exp(-d*TTM)/B;     % Forward in G&C Model

% Define method names for reference
methodNames={'by blkprice','CRR','MonteCarlo'}; 
prices=zeros(3,1); % Preallocate a vector to store the option prices computed by the three methods

% Loop over the 3 methods:
% 1) Closed-form formula (Black-Scholes, "blkprice")
% 2) Cox-Ross-Rubinstein (CRR) binomial tree
% 3) Monte Carlo

for i=1:3
    pricingMode = i; %
    M=100; % M = number of simulations for MC, steps for CRR;

    % Compute the European option price via a custom function EuropeanOptionPrice
    OptionPrice = EuropeanOptionPrice(F0,K,B,TTM,sigma,pricingMode,M,flag);
    prices(i)=OptionPrice;
end

% Create and display a table with the pricing results for each method
T=table(methodNames',prices,'VariableNames',{'Method', 'Price'});
disp(T); 

% Comments on errors for different methods:
% - For a European option in the CRR model, the error behaves like O(1/N). 
%   Increasing N reduces the error but increases computation time.
% - For Monte Carlo, the statistical error behaves like O(1/sqrt(M)). 
%   Increasing M reduces the estimator's standard deviation, again at higher computational cost.
% - You can estimate the standard error in a binomial or MC approach by computing:
%   payoffStd = std(payoffs);
%   stdError  = B * payoffStd / sqrt(M);

%% b) Errors Rescaling 

% Introduce the exact price via the closed-form Black-Scholes formula
C_exact = EuropeanOptionClosed(F0,K,B,TTM,sigma,flag);
fprintf('Exact price, closed formula = %.6f\n', C_exact);

% Evaluate the CRR error as the number of steps N increases
N_list = 2.^[1:10];  % We use powers of 2 for convenience
errorCRR = zeros(size(N_list));

bidAskThreshold = 1e-4;  % 1 bp = 0.0001


for i=1:length(N_list)
    N=N_list(i);

    % Compute the price via CRR with N steps
    C_CRR=EuropeanOptionCRR(F0,K,B,TTM,sigma,N,flag);

    % Calculate the absolute error compared to the exact BS price
    errorCRR(i)=abs(C_CRR-C_exact);

    % If the error is below the bid-ask threshold, we can stop increasing N
    if errorCRR(i)<=bidAskThreshold
        break
    end
end


% --- CRR Error Plot ---
figure('Name','Errore CRR vs N','NumberTitle','off');
semilogy(N_list, errorCRR, '-o','LineWidth',1.5);
hold on;
yline(bidAskThreshold, '--r', '1 bp','LineWidth',1.2);
xlabel('N (steps CRR)');
ylabel('Absolute error |CRR - Exact|');
title('Error in the CRR model');
grid on;
% From this plot, we can see at which N the error first drops below 1 bp.


%% Monte Carlo Standard Error Variation
% Evaluate how the Monte Carlo standard error decreases as M (the number of simulations) grows
M_list = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400,204800, 409600, 550000];
errorMC = zeros(size(M_list));

for i = 1:length(M_list)
    M = M_list(i);
    
    % Compute both the MC price and its standard error
    [C_MC, stdErrorMC] = EuropeanOptionMC_withStd(F0, K, B, TTM, sigma, M, flag);
    
    % we consider the unbiased standard deviation  
    % (i.e., the standard error of the MC estimator) as the "error"
    errorMC(i) = stdErrorMC;
end

% --- Monte Carlo Error Plot ---
figure('Name','Std Error MC vs M','NumberTitle','off');
semilogy(M_list, errorMC, '-o','LineWidth',1.5);
hold on;
yline(bidAskThreshold, '--r', '1 bp','LineWidth',1.2);
xlabel('M (MC simulations)');
ylabel('Std Error MC');
title('Error in the Monte Carlo model (Std Error)');
grid on;


%% c) Graphs in log-log scale for show how the numerical errors for a call rescale

% Plot in a log-log scale the CRR numerical errors wrt the number of steps M increases 
PlotErrorCRR(F0,K,B,TTM,sigma);

% Plot in a log-log scale the MC numerical errors wrt the number of simulations M increases
[M, stdEstim] = PlotErrorMC(F0,K,B,TTM,sigma);

%% d) KO Option
KO=1.4;
BarrierNames={'Closed Formula','CRR','MonteCarlo'}; %defining names of numerical methods
Barr=zeros(3,1);
%computing the prices both with numerical methods and closed formula 
BarrierPriceMC = EuropeanOptionKOMC(F0,K,KO,B,TTM,sigma,550000);
BarrierPriceClosed = EuropeanOptionKOClosed(F0,K,KO,B,TTM,sigma,flag);
BarrierPriceCRR = EuropeanOptionKOCRR(F0,K,KO,B,TTM,sigma,80);

%defining and plotting the table with obtained values
Barr=[BarrierPriceClosed; BarrierPriceCRR; BarrierPriceMC];
T2=table(BarrierNames', Barr,'VariableNames',{'Method', 'Price'});
disp(T2); 
%% e) KO Option Vega
%creating the vectors for plotting 
xAxis = linspace(0.65, 1.45, 100);
vegaVectorCRR= zeros(length(xAxis), 1);
vegaVectorMC= zeros(length(xAxis), 1);
vegaVectorExact = zeros(length(xAxis), 1);
%computing the values for plotting
for index = 1:length(xAxis)
    spot = xAxis(index);  % Use values from xAxes directly
    N=80;
    vegaVectorCRR(index) = vegaKO(spot/B, K, KO, B, TTM, sigma, N, 1);
    N = 500000;
    vegaVectorMC(index) = vegaKO(spot/B, K, KO, B, TTM, sigma, N, 2);
    vegaVectorExact(index) = vegaKO(spot/B, K, KO, B, TTM, sigma, N, 3);
end

% Plotting the results
figure();
plot(xAxis, vegaVectorCRR);
title('CRR Approximation of Vega');
xlabel('Spot Price');
ylabel('Vega');

figure();
plot(xAxis, vegaVectorMC);
title('MC Approximation of Vega');
xlabel('Spot Price');
ylabel('Vega');

figure();
plot(xAxis, vegaVectorExact);
title('Exact Vega');
xlabel('Spot Price');
ylabel('Vega');

%% f) Antithetic Variables
% Compare the MC method with the MC method using the antithetic variables technique

% Calculates the standard error of the MC estimate with AV
[M, stdEstimAV] = PlotErrorMCAntitheticVariables(F0,K,B,TTM,sigma);

%Compare plot of the errors
figure()
% Monte Carlo Error
loglog(M, stdEstim, 'o-','LineWidth',1.5,'DisplayName','MC StdError', Color='b');
hold on; grid on;
% Monte Carlo Error with AV
loglog(M, stdEstimAV, 'o-','LineWidth',1.5,'DisplayName','MC StdErrorAV',Color='r');
xlabel('M (simulation MC)');
ylabel('Std Error');
title('Error MC vs Error AV MC');
legend('Location','best');
hold off;

% Let's calculate the percentage of cases where stdEstimAV is better than stdEstim
% (stdEstimAV < stdEstim)
percentageAVbetter = (sum(stdEstimAV < stdEstim) / length(stdEstim)) * 100;

% Print result
fprintf('Percentage of cases where stdEstimAV is less than stdEstim: %.2f%%\n', percentageAVbetter);


%% g) Bermudan option CRR

N_opt = 80; % Number of steps in the binomial tree
flag = 1; % 1 for Call option, -1 for Put option

% Compute the price of the Bermudan option using the CRR method 
optionPriceBerm = BermudanOptionCRR(F0,K,B,TTM,sigma,N_opt,flag);

% Display the computed option price
fprintf('The Bermudian option price is %f\n', optionPriceBerm );



%% h) Bermudan and European option CRR varying the dividend yield

% Define a range of dividend yields from 0% to 5%
div=linspace(0,0.05,10);

optionPriceBerm_div=zeros(length(div),1);
optionPriceEU_div=zeros(length(div),1);

% Compute option prices for different dividend yields
for i=1:length(div)
    optionPriceBerm_div(i) = BermudanOptionCRR_div(F0,K,B,div(i),TTM,sigma,N_opt,flag);
    optionPriceEU_div(i) = EuropeanOptionCRR_div(F0,K,B,div(i),TTM,sigma,N_opt,flag);
end


% Plot: Option Prices vs Dividend Yield
figure();
plot(div * 100, optionPriceEU_div, 'b-o', 'LineWidth', 2); hold on;
plot(div * 100, optionPriceBerm_div, 'r-s', 'LineWidth', 2);
xlabel('Dividend Yield (%)');
ylabel('Option Price');
title('European vs. Bermudan Option Prices vs. Dividend Yield');
legend('European', 'Bermudan', 'Location', 'best');
grid on;

% Plot: Price Difference (Bermudan - European) vs Dividend Yield
figure();
plot(div * 100, optionPriceBerm_div - optionPriceEU_div, 'k-', 'LineWidth', 2);
xlabel('Dividend Yield (%)');
ylabel('Price Difference');
title('Price Premium: Bermudan - European');
grid on;


