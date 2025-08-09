function [datesCDS, survProbs, intensities] = bootstrapCDS_accrual(datesDF, discounts, datesCDS, spreadsCDS, recovery)
% bootstrapCDS_accrual computes the survival probabilities and hazard rates
% (intensities) for a set of CDS instruments, including an accrual term.
%
% INPUT:
%   datesDF   : Vector of dates corresponding to discount factors
%   discounts : Vector of discount factors
%   datesCDS  : Vector of CDS maturity/expiry dates
%   spreadsCDS: Vector of CDS spreads
%   recovery  : Recovery rate for the reference entity
%
% OUTPUT:
%   datesCDS   : Same as input, dates of expiry of the CDS instruments
%   survProbs  : Survival probabilities at each CDS date
%   intensities: Hazard rates (intensities) at each CDS date
%
% DESCRIPTION:
%   This function bootstraps the survival probabilities and hazard rates
%   using a piecewise-constant hazard rate model, explicitly including the
%   accrual on the premium leg of the CDS. The accrual term accounts for the
%   fraction of the premium that would be paid if default occurs between
%   payment dates.

%% --- PRELIMINARY COMPUTATIONS -------------------------------------------

% The settlement date is taken as the first date in the discount factor array.
settlement = datesDF(1);

% Convert discount factor dates to year fractions relative to 'settlement'
% using day-count convention #3 (e.g., Act/365 or similar).
dates_zeroRates = yearfrac(settlement, datesDF, 3);

% Compute zero rates from discount factors. 'zeroRates' is assumed to be 
% a function that extracts zero rates from discount factors.
zRates = zeroRates(datesDF, discounts);

% Convert CDS expiry dates to year fractions relative to 'settlement' 
% for interpolating zero rates.
dates_CDS_zrates = yearfrac(settlement, datesCDS, 3);

% Interpolate the zero rates at the CDS expiry dates and convert to
% discount factors. The factor "/100" assumes zRates is in percentage form.
discounts_CDS = exp( -interp1(dates_zeroRates, zRates, dates_CDS_zrates) ...
    / 100 .* dates_CDS_zrates );

% Compute time intervals (delta_t) for each CDS maturity date relative
% to the settlement, using day-count convention #6 (e.g., 30/360).
% The 'diff' of [0; yearfrac(...)] gives the incremental year fraction 
% between consecutive CDS maturities, with the first element referencing 0.
delta_t = diff([0; yearfrac(settlement, datesCDS, 6)]);

% Number of CDS maturities/spreads
N = length(spreadsCDS);

% Initialize output arrays for survival probabilities and intensities
survProbs = zeros(N,1);
intensities = zeros(N,1);

%% --- BOOTSTRAPPING PROCEDURE -------------------------------------------

% ------------------------------
% First iteration (i = 1) setup
% ------------------------------
%
% The formula for the first survival probability includes an accrual term 
% proportional to (delta_t(1)/2). The approach is an approximation where 
% the default is assumed to happen in the middle of the period on average.

survProbs(1) = (1 - recovery - spreadsCDS(1)*delta_t(1)/2) / ...
               (1 - recovery + spreadsCDS(1)*delta_t(1)/2);

% First hazard rate (intensity) from survival probability over delta_t(1)
intensities(1) = -log(survProbs(1)) / delta_t(1);

% Initialize quantities for the bootstrapping loop
% BPV (risky annuity) accumulates discounted survival probabilities
BPV = delta_t(1) * survProbs(1) * discounts_CDS(1);

% e accumulates the discounted expected losses (1 - survivalProb)
e = discounts_CDS(1) * (1 - survProbs(1));

% e_cum_deltas tracks the discounted accrual portion
e_cum_deltas = delta_t(1)/2 * discounts_CDS(1) * (1 - survProbs(1));

% --------------------------------------------------
% Iteration from i = 2 to N for subsequent maturities
% --------------------------------------------------
for i = 2:N
    
    % Numerator of the survival probability formula:
    %  (1 - recovery - s_i * delta_t(i)/2) * discount * S_{i-1}
    %  + (1 - recovery)*e - s_i*(e_cum_deltas + BPV)
    %
    % Denominator:
    %  (1 - recovery + s_i * delta_t(i)/2)*discount
    num = (1 - recovery - spreadsCDS(i)*delta_t(i)/2) * discounts_CDS(i) * survProbs(i-1) ...
        + (1 - recovery)*e ...
        - spreadsCDS(i)*e_cum_deltas ...
        - spreadsCDS(i)*BPV;
    
    den = (1 - recovery + spreadsCDS(i)*delta_t(i)/2) * discounts_CDS(i);
    
    % Compute the survival probability for maturity i
    survProbs(i) = num / den;
    
    % Update BPV: add the present value of the new 'slice' of survival
    BPV = BPV + delta_t(i) * survProbs(i) * discounts_CDS(i);
    
    % Update e: add the difference in survival probability discounted
    e = e + discounts_CDS(i) * (survProbs(i-1) - survProbs(i));
    
    % Update e_cum_deltas: accrual part
    e_cum_deltas = e_cum_deltas + (delta_t(i)/2) * discounts_CDS(i) * ...
        (survProbs(i-1) - survProbs(i));
    
    % Compute hazard rate for the interval [i-1, i]
    % The formula uses log of the ratio of consecutive survival probabilities.
    intensities(i) = -log( survProbs(i) / survProbs(i-1) ) / delta_t(i);
    
end

end
