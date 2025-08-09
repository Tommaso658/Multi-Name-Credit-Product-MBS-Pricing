function npv_value = computeNPV3(dates, discounts,settelment_date,first_cashflow_date, C0, AAGR, NA)
    % computeNPV - Computes the NPV of monthly cash flows that grow annual
    % dates               = dates obtained by boostrap
    % discounts           = dates obtained by boostrap  
    % settelment_date     = settelment date
    % first_cashflow_date = first cash_flow date
    % C0                  = initial cash flow (in euros)
    % AAGR                = annual growth rate (e.g., 0.05 for 5%)
    % NA                  = total number of years

    %Vector with all date of cashflow payment
    date_25 = busdate(datenum(first_cashflow_date:calmonths(1):first_cashflow_date + calyears(20))-1);
    
    %Interpolation for found discount factor
    zeroRates_discount = interp1(dates, zeroRates(dates, discounts), date_25, 'linear');
    %Calculate discount factor da zeroRates
    discount_factors = exp(-yearfrac(settelment_date, date_25(1:end),3).*(zeroRates_discount/100));

    % Month index (from 0 to 240)
    t = (0:NA*12);
    
    % Number of years passed for each month (k = floor(t / 12))
    k = floor(t / 12);
    
    % Vectorized computation of cash flows with annual growth
    cash_flows = C0 * (1 + AAGR) .^ k;
    
    % Vectorized NPV computation
    npv_value = sum(cash_flows .* discount_factors);

end