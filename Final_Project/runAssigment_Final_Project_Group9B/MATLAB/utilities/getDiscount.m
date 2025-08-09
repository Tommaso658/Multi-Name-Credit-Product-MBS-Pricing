function new_discounts = getDiscount(discounts,new_dates,dates,settlement_date)
    act365 = 6;
    rates = zeroRates(dates,discounts);
    new_rates = interp1(dates,rates,new_dates);
    new_discounts = exp(-new_rates/100.*yearfrac(settlement_date,new_dates,act365));
end