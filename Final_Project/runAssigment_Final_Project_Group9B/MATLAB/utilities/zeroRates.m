function zRates = zeroRates(dates, discounts)

% computation of zero rates knowing discount factors
% y(t) = - ln(B(t0,t)) / delta(t0,t) (with delta in Act365)

act365 = 6;
delta = yearfrac(dates(1),dates(2:end),act365);
zRates = [0; -log(discounts(2:end))./delta.*100]; %percentage

end