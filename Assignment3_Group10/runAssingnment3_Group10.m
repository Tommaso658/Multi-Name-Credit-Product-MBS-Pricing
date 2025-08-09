% runAssignment3
% Group 10, Academic Year 2024-2025
% This script executes the assignment, covering market data bootstrapping,
% asset swap spread computation, CDS bootstrapping, credit simulation, and
% MBS pricing.

clear all;      % Clear workspace variables
close all;      % Close all open figures
clc;            % Clear command window

%% Settings
formatData = 'dd/mm/yyyy';  % Date format string – adjust based on computer settings

%% Read Market Data
% Reads market data from an Excel file. Note: this function is designed
% for Windows OS; modifications may be required for other operating systems.
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap', formatData);

%% Bootstrap the Discount Curve
% The bootstrap function returns a set of dates (including settlement) and
% the corresponding discount factors computed from the market data.
[dates, discounts] = bootstrap(datesSet, ratesSet);
% Compute zero rates from the obtained discount factors.
z_Rates = zeroRates(dates, discounts);

%% 1) Asset Swap Spread Euribor3m

% Define input data for the asset swap
FV = 1;                      % Face value (normalized to 1)
clean_price = 1.01;          % Clean price (101% of face value)
expiry = 5;                  % Maturity in years
coupon = 4.8/100;            % Annual coupon rate (4.8%)
issue_date = datenum('31-Mar-2022');      % Issue date
settlement_date = datenum('02-Feb-2023');   % Settlement date

% Compute the accrued interest from issue to settlement using day count
% convention 6 (e.g., 30/360).
accrual = coupon * yearfrac(issue_date, settlement_date, 6);
% Dirty price is the sum of the clean price and the accrued interest.
dirty_price = clean_price + accrual;

% Define payment dates for the fixed leg (annual payments for 6 years)
% adjust for business days:
payament_dates_fixed_leg = ['31-Mar-2023';'31-Mar-2024';'31-Mar-2025';'31-Mar-2026';'31-Mar-2027';'31-Mar-2028'];
payament_dates_fixed_leg = datenum(busdate(payament_dates_fixed_leg),1);

% Define payment dates for the floating leg:
% The first payment occurs on 31-Mar-2023 and subsequent payments are every 3 months.
first_float_payament = datenum('31-Mar-2023');
num_float_payments = expiry * 12/3 + 1;  % Total number of floating payments
payament_dates_floating_leg = zeros(num_float_payments, 1);
payament_dates_floating_leg(1) = first_float_payament;

% Loop to compute the remaining floating payment dates
for ii = 2:num_float_payments - 1
    % Add 3 months to the previous payment date
    payament_dates_floating_leg(ii) = addtodate(payament_dates_floating_leg(ii-1), 3, "month");
    % Adjust payment date if it falls on a weekend (Saturday=7 or Sunday=1)
    Day = weekday(payament_dates_floating_leg(ii));
    if Day == 1 || Day == 7
        payament_dates_floating_leg(ii) = busdate(payament_dates_floating_leg(ii), 1);
    end
end
% Ensure the last floating payment date is set to 31-Mar-2028
payament_dates_floating_leg(end) = datenum('31-Mar-2028');

% Compute year fractions from settlement date for discounting purposes
yf_dates = yearfrac(settlement_date, dates, 3);

% Compute discount factors for the floating leg payments
dates_discounts_float = yearfrac(settlement_date, payament_dates_floating_leg, 3);
discounts_factor_float = exp((-interp1(yf_dates, z_Rates, dates_discounts_float)/100) .* dates_discounts_float);

% Compute discount factors for the fixed leg payments
dates_discounts_fixed = yearfrac(settlement_date, payament_dates_fixed_leg, 3);
discounts_factor_fixed = exp((-interp1(yf_dates, z_Rates, dates_discounts_fixed)/100) .* dates_discounts_fixed);

% Calculate time increments (delta_t) for basis point value (BPV) calculation
delta_t_float = diff([0; yearfrac(settlement_date, payament_dates_floating_leg, 6)]);
delta_t_fixed = diff([0; yearfrac(settlement_date, payament_dates_fixed_leg, 6)]);

% Compute Basis Point Value (BPV) for fixed and floating legs
BPV_fixed = sum(delta_t_fixed .* discounts_factor_fixed);
BPV_float = sum(delta_t_float .* discounts_factor_float);

% Compute the price of the fixed leg: last discount factor plus coupon payments
price = discounts_factor_fixed(end) + coupon * BPV_fixed;

% Compute the asset swap spread (s_asw) in terms of BPV
s_asw = (price - dirty_price) / BPV_float;
fprintf("The s_asw is: %.6f bps \n\n", s_asw * 10000); % Multiply by 10,000 to convert to basis points

%% 2) CDS Bootstrap

% Set recovery rate for the CDS
recovery = 40/100;
% Define CDS spreads (in decimal) for maturities corresponding to datesSet.swaps
spreadsCDS = [40 44 47 49 51 52]' * 1e-4;
% Extract the CDS maturity dates from the market data structure (first 6 entries)
datesCDS = datesSet.swaps(1:6);

% Compute CDS bootstrapping using different methods:

% 1) Approximation neglecting the accrual term (flag = 1)
flag = 1;
[datesCDS, survProbs_Noaccural, intensities_Noaccrual] = bootstrapCDS(dates, discounts, datesCDS, spreadsCDS, flag, recovery);

% 2) Exact method including the accrual term (flag = 2)
flag = 2;
[datesCDS, survProbs_accrual, intensities_accrual] = bootstrapCDS(dates, discounts, datesCDS, spreadsCDS, flag, recovery);

% 3) Jarrow–Turnbull (JT) approach (flag = 3)
flag = 3;
[datesCDS, survProbs_JT, intensities_JT] = bootstrapCDS(dates, discounts, datesCDS, spreadsCDS, flag, recovery);

% Display the maximum differences between the methods for survival probabilities and intensities
fprintf("Max prob error neglecting accrual term: %.2d \n", max(abs(survProbs_accrual - survProbs_Noaccural)) );
fprintf("Max intensities error neglecting accrual term: %.2d \n", max(abs(intensities_accrual - intensities_Noaccrual)) );
fprintf("\nMax prob error with JT approximation: %.2d \n", max(abs(survProbs_accrual - survProbs_JT)) );
fprintf("Max intensities error with JT approximation: %.2d \n\n", max(abs(intensities_accrual - intensities_JT)) );

% Plot survival probabilities for the three methods
figure;
x = (0:6);  % Time index: 0 represents time 0 (settlement), then each year
plot(x, [1; survProbs_Noaccural], '-s', 'LineWidth', 1.5, 'MarkerSize', 10);
hold on;
plot(x, [1; survProbs_accrual], '-o', 'LineWidth', 1.1, 'MarkerSize', 5);
plot(x, [1; survProbs_JT], '-^', 'LineWidth', 1.1, 'MarkerSize', 5);
legend('Approximate', 'Exact', 'JT', 'Location', 'northeast');
set(gca, 'YGrid', 'on');
title('Survival Probabilities', 'FontSize', 14);
xlabel('Time', 'FontSize', 12);
ylabel('Probability', 'FontSize', 12);
xticklabels({'02/02/23','1y','2y','3y','4y','5y','6y'});
xtickangle(45);
hold off;

% Plot intensities for the approximate and exact methods
figure;
x_i = linspace(0, 6, 500);
% Interpolate intensities for a smoother curve; using "next" interpolation method
plot( x_i, interp1(x, [intensities_Noaccrual(1); intensities_Noaccrual], x_i, "next", "extrap"), '-s', 'LineWidth', 1.5, 'MarkerSize', 0.75);
hold on;
plot( x_i, interp1(x, [intensities_accrual(1); intensities_accrual], x_i, "next", "extrap"), '-o', 'LineWidth', 1.1, 'MarkerSize', 0.55);
plot(x_i,interp1(x,[intensities_JT(1);intensities_JT],x_i,"next","extrap"),...
    '-^', 'LineWidth', 1.1, 'MarkerSize', 0.55);
legend('Approximate', 'Exact','JT ','Location', 'southeast');
set(gca, 'YGrid', 'on');
title('Intensities', 'FontSize', 14);
xlabel('Time', 'FontSize', 12);
ylabel('Intensity', 'FontSize', 12);
xticklabels({'02/02/23','1y','2y','3y','4y','5y','6y'});
xtickangle(45);
hold off;

%% 3) Credit Simulation

% Fix the random seed for reproducibility
rng(0)
% Define simulation parameters:
theta = 4;         % Threshold time separating two regimes
lambda1 = 5e-4;    % Hazard rate (intensity) before theta
lambda2 = 9e-4;    % Hazard rate (intensity) after theta
t0 = 0;            % Start time
T = 30;            % Final time horizon (30 years)
M = 1e5;           % Number of simulations

% Generate M uniform random numbers for simulation
u = rand(1, M);

% Simulate default times (tau) using a piecewise model:
% For u >= exp(-lambda1*(theta-t0)), use the first regime,
% otherwise use the second regime.
tau = (-log(u) / lambda1 + t0) .* (u >= exp(-lambda1*(theta-t0))) + ...
      (theta*(lambda2 - lambda1) + lambda1*t0 - log(u)) / lambda2 .* (u < exp(-lambda1*(theta-t0)));

% Plot simulated default times (tau) vs. uniform random numbers (u) on a semilog scale
figure();
semilogy(tau, u, 'o');
title('Simulated Survival Probability');
xlabel('tau');
ylabel('P(0,tau)');

% Count number of defaults in each regime for MLE estimation
n1 = sum(tau <= theta);   % Number of defaults before theta
n2 = M - n1;              % Number of defaults after theta

% Maximum Likelihood Estimators (MLE) for lambda1 and lambda2
L1_mle = n1 / (sum(tau .* (tau <= theta)) + n2 * theta);
L2_mle = n2 / (sum((tau - theta) .* (tau > theta)));

% Compute the information matrix elements for the MLEs
I = zeros(2);
I(1,1) = (1 - exp(-L1_mle * theta)) / L1_mle^2;
I(2,2) = exp(-L1_mle * theta) / L2_mle^2;

a = 0.95; % Confidence level (alpha)
% Confidence Intervals (CI) for the MLEs using the standard normal quantile
CI1 = [L1_mle - norminv((1 + a) / 2) / sqrt(M * I(1,1)), L1_mle + norminv((1 + a) / 2) / sqrt(M * I(1,1))];
CI2 = [L2_mle - norminv((1 + a) / 2) / sqrt(M * I(2,2)), L2_mle + norminv((1 + a) / 2) / sqrt(M * I(2,2))];

% Compute the lengths of the confidence intervals
length_CI1 = abs(diff(CI1));
length_CI2 = abs(diff(CI2));

% Print the confidence intervals (scaled by 10^4 for basis points)
fprintf("CI_lambda1: [%.4f, %.4f] x 10^-4\n", CI1 * 1e4);
fprintf("CI_lambda2: [%.4f, %.4f] x 10^-4\n", CI2 * 1e4);

%% Alternative Empirical Method for Survival Probability Estimation

% Sort the simulated default times
tau_s = sort(tau);
% Compute empirical survival probabilities by counting defaults
prob = (cumsum(ones(1, M), 'reverse') - ones(1, M)) / M;
% Separate the probabilities before and after the threshold theta
prob1 = prob(tau_s <= theta);
dates1 = tau_s(tau_s <= theta);  % Default times up to theta

% For times between theta and T, extract nonzero probabilities and corresponding times
prob2 = nonzeros(prob(tau_s > theta & tau_s <= T))';
dates2 = nonzeros(tau_s(tau_s > theta & tau_s <= T))';

% Plot empirical survival probabilities on a semilog scale
figure(3);
semilogy(dates1, prob1, '*', dates2, prob2, '*');
legend('Survival Probability before \theta', 'Survival Probability after \theta');
title('Comparison (Empirical Method)');
xlabel('\tau');
ylabel('P(0,\tau)');

% Plot simulated survival probabilities up to 30 years
figure(4);
semilogy(tau(tau <= T), u(tau <= T), '*');
title('Simulated Survival Probability up to 30 years');
xlabel('\tau');
ylabel('P(0,\tau)');

%% 4) MBS Pricing

% Pricing of a mezzanine tranche for ABS using three methods:
% HP (exact integration), KL (Kullback-Leibler based integration), and LHP.

% Define required variables for MBS pricing:
p = 0.05;            % Default probability for each mortgage
rho = 0.40;          % Correlation between mortgage defaults
pi = 0.20;           % Average recovery rate
Kd = 0.05;           % Lower detachment point: tranche starts losing when losses exceed this percentage
Ku = 0.09;           % Upper detachment point: tranche is completely wiped out beyond this loss percentage
I = linspace(10, 2e4, 2e4-9);  % Vector representing the number of obligors in the reference portfolio
notional = 1e9;      % Notional value of the portfolio
T = 3;               % Maturity in years

Nt = notional * (Ku - Kd);  % Notional of the tranche

discount = 0.914627080074705;  % 3-year discount factor (obtained from market data)

% Define sample values for the number of mortgages for HP and KL integration methods
I_HP = [10, 20, 40, 70, 100, 400, 700, 1000];   % Values of I for HP method
I_KL = [10, 20, 40, 70, 100, 400, 700, 1000, 4000, 7000, 10000, 20000];   % Values of I for KL method

Width = 1.5;   % Line width for plots

%% Mezzanine Tranche Pricing

% Exact Integration Method (HP)
pricesExact = zeros(length(I_HP), 1);
index = 1;
for i = I_HP
    % Loss function computes the expected loss for the given portfolio size
    pricesExact(index, 1) = 1 - lossFunction(i, pi, Ku, Kd, p, rho);
    index = index + 1;
end
% Discount the computed tranche prices
pricesExact = discount * pricesExact;

% Kullback-Leibler (KL) Integration Method
pricesKL = zeros(length(I_KL), 1);
index = 1;
for i = I_KL
    pricesKL(index, 1) = 1 - lossMKL(i, pi, Ku, Kd, p, rho);
    index = index + 1;
end 
pricesKL = discount * pricesKL;

% LHP Method (a simplified closed-form pricing)
lhp = ones(length(I_KL), 1) * LHP(pi, Ku, Kd, p, rho) * discount;

% Plot the mezzanine tranche price comparison using the three methods
figure(1)
semilogx(I_HP, pricesExact, 'x-', 'LineWidth', Width);
hold on
grid on
semilogx(I_KL, pricesKL, 'o-', 'LineWidth', Width);
semilogx(I_KL, lhp, 'LineWidth', Width);
title("Mezzanine Tranche Price Comparison");
xlabel("Number of Mortgages");
ylabel("% Price Tranche");
legend("Price Exact", "Price KL", "Price LHP", 'Location', 'best', 'FontSize', 12);
hold off

%%  Price EQUITY Tranche

% Adjust detachment points for the equity tranche:
Kd = 0;     % Equity tranche starts losing immediately
Ku = 0.05;  % Equity tranche is wiped out when losses exceed 5%

% Compute equity tranche prices using the Exact Integration Method (HP)
pricesExact = zeros(length(I_HP), 1);
index = 1;
for i = I_HP
    pricesExact(index, 1) = 1 - lossFunction(i, pi, Ku, Kd, p, rho);
    index = index + 1;
end
pricesExact = discount * pricesExact;

% Compute equity tranche prices using the KL Integration Method
pricesKL = zeros(length(I_KL), 1);
index = 1;
for i = I_KL
    pricesKL(index, 1) = 1 - lossMKL(i, pi, Ku, Kd, p, rho);
    index = index + 1;
end 
pricesKL = discount * pricesKL;

% LHP pricing for equity tranche (same approach as before)
lhp = ones(length(I_KL), 1) * LHP(pi, Ku, Kd, p, rho) * discount;

% Plot the equity tranche price comparison
figure(2)
semilogx(I_HP, pricesExact, 'x-', 'LineWidth', Width);
hold on
grid on
semilogx(I_KL, pricesKL, 'o-', 'LineWidth', Width);
semilogx(I_KL, lhp, 'LineWidth', Width);
title("Equity Tranche Price Comparison");
xlabel("Number of Mortgages");
ylabel("% Price Tranche");
legend("Price Exact", "Price KL", "Price LHP", 'Location', 'best', 'FontSize', 12);
hold off

%% EQUITY - Correction for KL Integration
% For the equity tranche, adjust the KL pricing method by considering the
% loss of the reference portfolio.

LossReferencePortfolio = (1 - pi) * p;  % Expected loss of the reference portfolio

% For the mezzanine tranche pricing using the corrected KL method, adjust the detachment points:
Kd = 0.05;   % Mezzanine tranche starts losing when losses exceed 5%
Ku = 1.00;   % Mezzanine tranche is completely wiped out when losses reach 100%
Nt = notional * (Ku - Kd);   % Notional of the mezzanine tranche

pricesKLEquity = zeros(length(I_KL), 1);
index = 1;
lossRP = zeros(length(I_KL), 1);
expectedLoss = zeros(length(I_KL), 1);

for i = I_KL
    % Compute the loss on the reference portfolio for a given number of obligors
    lossRP(index, 1) = Nt * lossMKL(i, pi, Ku, Kd, p, rho);
    % Expected loss on the entire portfolio
    expectedLoss = notional * (1 - pi) * p;
    % Adjusted price using the corrected KL method
    pricesKLEquity(index, 1) = (1 - (expectedLoss - lossRP(index, 1)) / (notional - Nt)) * discount;
    index = index + 1;
end

% Plot the equity tranche price comparison with the corrected KL method
figure(3)
semilogx(I_HP, pricesExact, 'x-', 'LineWidth', Width);
hold on
grid on
semilogx(I_KL, pricesKLEquity, 'o-', 'LineWidth', Width);
semilogx(I_KL, lhp, 'LineWidth', Width);
title("Equity Tranche Price Comparison (Corrected KL)");
xlabel("Number of Mortgages");
ylabel("% Price Tranche");
legend("Price Exact", "Price KL (Corrected)", "Price LHP", 'Location', 'best', 'FontSize', 12);
hold off
