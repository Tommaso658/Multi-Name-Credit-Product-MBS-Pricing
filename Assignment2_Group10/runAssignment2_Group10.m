% runAssignment2
%  group 10, AY2024-2025
% Computes Euribor 3m bootstrap with a single-curve model
%
% This is just code structure: it should be completed & modified (TBM)
%
% to run:
% > runAssignment2_TBM

clear all;
close all;
clc;
warning off;

%% Settings
formatData='dd/mm/yyyy'; %Pay attention to your computer settings 

%% Read market data
% This fuction works on Windows OS. Pay attention on other OS.

[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap', formatData);

%% 1) Bootstrap
% dates includes SettlementDate as first date

[dates, discounts] = bootstrap(datesSet, ratesSet);


%% Compute Zero Rates
z_Rates = zeroRates(dates, discounts);

%% Plot Results

figure(1);
title('reconstruction of the curves')
yyaxis left
plot(dates, discounts,'-o','LineWidth',1.5, 'MarkerSize',2);
xlabel('dates'); ylabel('discount factor');
hold on
yyaxis right
plot(dates, z_Rates*100,'*-','LineWidth',1.5,'MarkerSize',4);
ylabel('zero rate %');
datetick
set(gca, 'YGrid', 'on','XGrid','off')

%% 2) Pricing a 7y I.B. coupon bond with coupon rate = mid-mkt 7y swap rate

%the fixed coupon rate is the mid-mkt rate of the 7y swap
coupon_rate = 1/2*(ratesSet.swaps(6,1) + ratesSet.swaps(6,2));
face_value = 10^8; %FV = 100Mln
settlement = datesSet.settlement;

%the first coupon payament happens one year after the settlement date
year_one = settlement + 365;

%we check if the first payament is a buisness day
%if not we go one step ahead
while isweekend(datetime(datestr(year_one)))  
    year_one = year_one + 1; 
end

%dates in which the coupon payaments happen
coupon_dates = [year_one; datesSet.swaps(1:6)];

i = find(dates > coupon_dates(1), 1);
zrates = zeroRates([settlement dates(i-1) dates(i)], [1 discounts(i-1) discounts(i)]);
zrates = zrates/100;

%zero rates interpolation in order to get the discount factor at 1y
zrates_1y = interp1([dates(i-1) dates(i)], zrates(2:end), coupon_dates(1),'spline'); 
 
discount_1y = exp(-zrates_1y * yearfrac(settlement,coupon_dates(1),3)); %from zero rates to discount factor at 1y
discounts_coupon = [discount_1y discounts(13:18)];    


%year fraction with convention 30/360E
year_frac = yearfrac([settlement; coupon_dates(1:end-1)],coupon_dates(1:end),6);

cash_flows = coupon_rate * year_frac * face_value;
cash_flows(end) = cash_flows(end) + face_value; % at maturity we have to include also the face value



% price as the actualized sum of the cash flows
price_7y_coupon_bond = sum(cash_flows .* discounts_coupon');

fprintf('The price of the 7y_coupon_bond is %e\n', price_7y_coupon_bond);

%% 3) DV01 for an IRS, Modified duration for a coupon bond

% Parameters swap 7 years
notional = 1e8; 
FixedRate = 0.028175;  % 2.8175% 
setDate = datesSet.settlement; %settlement date
shift_bps = 1e-4; % 1 basis point = 0.01% 

% Create the shifted rates from the initial market rates
ratesSetshift = ratesSet;
ratesSetshift.depos = ratesSetshift.depos + shift_bps;
ratesSetshift.futures = ratesSetshift.futures + shift_bps;
ratesSetshift.swaps = ratesSetshift.swaps + shift_bps;

%creation of the new curve with rates shifted
[dates, discounts_DV01] = bootstrap(datesSet, ratesSetshift);
dates = dates';
discounts_DV01 = discounts_DV01';
discounts = discounts';


one_year = setDate + 365;
while isweekend(datetime(datestr(one_year)))   
    one_year = one_year +1;     %check whether it is a weekend          
end
fixedLegDates = [one_year; datesSet.swaps(1:6)]; %dates of fixed leg, one each year starting from settlement date, for 7y


[DV01, BPV, DV01_z] = sensSwap( setDate, ...
                                fixedLegDates, ...
                                FixedRate, ...
                                dates, ...
                                discounts, ...
                                discounts_DV01);

notionals = [DV01 BPV DV01_z] * notional;


%Macaulay Duration 
MacD = sensCouponBond( setDate, fixedLegDates, FixedRate, dates', discounts');


% print the results
fprintf('--- results: ---\n');
fprintf('DV01 Swap 7y     = %.4e \n', DV01);
fprintf('DV01_Z Swap 7y     = %.4e \n',  DV01_z);
fprintf('BPV  Swap 7y     = %.4e \n',  BPV);
fprintf('Macaulay Dur Bond = %.4f years\n', MacD);

%% 4) Net Present Value (NPV)

%Initial data
settelment_date = datetime('02/02/2023','InputFormat','dd/MM/yyyy');   %Settelment Date
first_cashflow_date = datetime('25/09/2026','InputFormat','dd/MM/yyyy');   %First cash_flow date
N = 20;  % Total number of cash flow months (20 years = 240 months)
AAGR = 0.05;   %Annual growth rate
C0 = [1500, 6000];   %Initial CashFlow


%NPV value from C0=1500
npv_value1 = computeNPV3(dates, discounts,settelment_date,first_cashflow_date, C0(1), AAGR, N)
%NPV value from C0=3000
npv_value2 = computeNPV3(dates, discounts,settelment_date,first_cashflow_date, C0(2), AAGR, N)
%Ratio between the two NPV
Ratio = npv_value2/npv_value1











