function [dates, discounts]=bootstrap(datesSet, ratesSet)
% Performs bootstrapping to construct a discount curve from market instruments
%
% INPUT:
% datesSet:  structure containing the settlement date and maturities of market instruments
%            - datesSet.settlement: settlement date
%            - datesSet.depos: maturities of deposit rates
%            - datesSet.futures: maturities of STIR futures
%            - datesSet.swaps: maturities of swap rates
%
% ratesSet:  structure containing bid and ask rates for market instruments
%            - ratesSet.depos: bid-ask rates for deposits
%            - ratesSet.futures: bid-ask rates for futures
%            - ratesSet.swaps: bid-ask rates for swaps
%
% OUTPUT:
% dates:      vector of maturity dates for the bootstrapped discount factors
% discounts:  vector of bootstrapped discount factors B(t0, ti)
%
% DESCRIPTION:
% This function bootstraps a discount curve using deposit rates, STIR futures, 
% and swap rates. The process follows these steps:
% 1. Adds the settlement date as the initial point in the curve.
% 2. Computes discount factors for deposit instruments.
% 3. Computes discount factors for STIR futures using interpolation/extrapolation.
% 4. Computes discount factors for swap instruments using iterative calculation.
%%
format long
settlement = datesSet.settlement;
dates = settlement; % we start setting as first date the settlement
discounts = 1; % first discout factor B(t0,t0) 
act360 = 2;
thirty360 = 6;
act365=3;
settlement=change_convention_date(settlement);

%% Deposits
first_fut_settle = change_convention_date(datesSet.futures(1,1)); % settle of the 1st future
dates_after_first_settle = find(change_convention_date(datesSet.depos) > first_fut_settle); % the position of the first depo to not consider 
n_depos = [1:dates_after_first_settle(1)]; % number of depos to consider
dates = [dates; change_convention_date(datesSet.depos(n_depos))]; % add depos' dates to the settlement

mid_depos = (ratesSet.depos(n_depos,1) + ratesSet.depos(n_depos,2))/2; % MID = (BID+ASK)/2
delta_depos = yearfrac(settlement,datesSet.depos(n_depos),act360); % delta(t0,ti)
disc_depos = 1./(1+delta_depos.*mid_depos); % B(t0,ti)

discounts = [discounts; disc_depos];

%% STIR futures
% B(t0, ti2) = B(t0, ti1) * B(t0; ti2, ti2)
% B(t0; ti1, ti2) = 1 / (1 + delta(ti, ti1)*L(t0; ti1, ti2)) 
% L(t0; ti1, ti2) is observed on the market

% using the values we know
mid_futures = (ratesSet.futures(:,1)+ratesSet.futures(:,2))/2; % MID = (BID+ASK)/2
delta_fut = yearfrac(datesSet.futures(:,1),datesSet.futures(:,2),act360); % delta(t0,ti)
fwd_discounts = 1./(1+delta_fut.*mid_futures); % B(t0; ti1, ti2))

% we start interpolation
n_fut = 8; % we take the most liquid (data from the excel spreadsheet)

y0= -log(discounts(end))/yearfrac(settlement,dates(end),act365);
delta1= yearfrac(settlement,datesSet.depos(length(n_depos)+1),act360);
mid_dep1= (ratesSet.depos(length(n_depos)+1,1) + ratesSet.depos(length(n_depos)+1,2))/2;
disc1=1/(1+delta1.*mid_dep1);
y1=-log(disc1)/yearfrac(settlement,datesSet.depos(length(n_depos)+1),act365);
y= y0 + (y1-y0)/(yearfrac(dates(end),datesSet.depos(length(n_depos)+1)))*(yearfrac(dates(end),datesSet.futures(1,1)));
B_interp=exp(-y*yearfrac(settlement,datesSet.futures(1,1)));
B_expiry = B_interp * fwd_discounts(1);
discounts=[discounts;B_expiry];
dates=[dates;change_convention_date(datesSet.futures(1,2))];
zero_rates= zeroRates(dates,discounts);
z_sett=interp1(dates,zero_rates,datesSet.futures(2,1),"linear","extrap")/100;
B_sett=exp(-z_sett*yearfrac(settlement,datesSet.futures(2,1)));
dates=[dates;change_convention_date(datesSet.futures(2,2))];
B_expiry_2=B_sett*fwd_discounts(2);
discounts=[discounts;B_expiry_2];

for i=3:n_fut
    
    settle = change_convention_date(datesSet.futures(i,1)); % ti1
    expiry = datesSet.futures(i,2); % ti2
    
    % interpolation/extrapolation
    B_settle = discounts(end);  % B(t0, ti1)
    
    % B(t0, ti2) = B(t0, ti1) * B(t0; ti1, ti2) 
    B_expiry = B_settle * fwd_discounts(i);
    
    dates = [dates; expiry];
    discounts = [discounts; B_expiry];

end

%% Swaps
% initialization
y0 = -log(discounts(end-1))/yearfrac(settlement,dates(end-1),act365);
y1 = -log(discounts(end))/yearfrac(settlement,dates(end),act365);
y = y0 + (y1-y0)*(datesSet.swaps(1)-dates(end-1))/(dates(end)-dates(end-1));
B = exp(-y*yearfrac(datesSet.settlement,datesSet.swaps(1),act365));
B01 = interp1(dates,discounts,datesSet.swaps(1),"linear"); % B(t0,t1), first year done with interp
BPV = B01*yearfrac(settlement, datesSet.swaps(1),thirty360); % initial BPV = B(t0,t1)*delta(t0,t1)
mid_swaps = (ratesSet.swaps(:,1)+ratesSet.swaps(:,2))/2; % MID = (BID+ASK)/2
[swap_date,swap_rate]=NewDates(datesSet.swaps,mid_swaps);
last_futures_expiry = datesSet.futures(n_fut,2); % expiry of the last future considered
first_swap = find(datesSet.swaps > last_futures_expiry, 1); % 1 stands for "first swap to consider"

n_swaps = length(swap_date);
for i=first_swap:n_swaps
    delta_0i = yearfrac(swap_date(i-1), swap_date(i), thirty360); % delta(t0;ti-1;ti)
    B0i = (1 - swap_rate(i)*BPV)/(1 + delta_0i*swap_rate(i)); % B(t0,ti); midswaps(i))=S0
    BPV = BPV + B0i*delta_0i; % BPV added with the new discount
    
    dates = [dates; swap_date(i)];
    discounts = [discounts; B0i];
    
end
end