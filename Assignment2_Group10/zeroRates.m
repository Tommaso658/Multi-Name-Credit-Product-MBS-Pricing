function zRates = zeroRates(dates, discounts)
%this function gives in output the zero rates at time indicated by dates

% INPUT:
% dates:        set of the dates
% discounts:    discount factor B(0,dates) 

% OUTPUT:
% zRates:       zero rates at time ti in percentege unit


% dates convertion for zRates = Act/365
zRates_conv = 3;
yf = yearfrac(dates(1), dates, zRates_conv);
zRates = -log(discounts)./yf;
zRates(1) = 0; % the zrate in t0 = 0
zRates = zRates * 100; %percentege unit

end