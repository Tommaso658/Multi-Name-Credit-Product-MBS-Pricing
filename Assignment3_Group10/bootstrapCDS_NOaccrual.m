function [datesCDS, survProbs, intensities] = bootstrapCDS_NOaccrual(datesDF, discounts, datesCDS, spreadsCDS, recovery)
% bootstrapCDS_NOaccrual computes survival probabilities and hazard rates 
% (intensities) for a set of CDS contracts, neglecting the accrual on the premium leg.
%
% INPUT:
%   datesDF    : Vector of dates corresponding to discount factors.
%   discounts  : Vector of discount factors.
%   datesCDS   : Vector of CDS expiry/maturity dates.
%   spreadsCDS : Vector of CDS spreads.
%   recovery   : Recovery rate for the CDS.
%
% OUTPUT:
%   datesCDS   : The same CDS expiry dates as input.
%   survProbs  : Computed survival probabilities at each CDS date.
%   intensities: Hazard rates (intensities) corresponding to each CDS date.

%% --- PRELIMINARY CALCULATIONS -------------------------------------------
% Define the settlement date as the first discount factor date.
settlement = datesDF(1);

% Convert discount factor dates to year fractions from settlement using
% day count convention 3 (e.g., Act/365 or similar).
dates_zeroRates = yearfrac(settlement, datesDF, 3);

% Compute zero rates from the discount factors.
zRates = zeroRates(datesDF, discounts);

% Convert CDS expiry dates into year fractions relative to the settlement date.
dates_CDS_zrates = yearfrac(settlement, datesCDS, 3);

% Interpolate zero rates at the CDS expiry dates and convert them into 
% discount factors. The division by 100 assumes the zero rates are expressed in percentage.
discounts_CDS = exp(-interp1(dates_zeroRates, zRates, dates_CDS_zrates)/100 .* dates_CDS_zrates);

% Compute the time increments (delta_t) for each CDS maturity.
% The 'diff' of the array [0; yearfrac(...)] gives the time length of each period 
% using day count convention 6 (e.g., 30/360).
delta_t = diff([0; yearfrac(settlement, datesCDS, 6)]);

%% --- INITIALIZATION -------------------------------------------------------
% Determine the number of CDS maturities/spreads.
N = length(spreadsCDS);

% Preallocate arrays for survival probabilities and hazard intensities.
survProbs = zeros(N, 1);
intensities = zeros(N, 1);

%% --- BOOTSTRAP PROCEDURE --------------------------------------------------
% First iteration (i = 1): Compute the initial survival probability.
% The formula used neglects the accrual term and adjusts only for recovery and spread.
survProbs(1) = (1 - recovery) / (1 - recovery + spreadsCDS(1) * delta_t(1));

% Calculate the first hazard rate (intensity) from the initial survival probability.
intensities(1) = -log(survProbs(1)) / delta_t(1);

% Initialize the cumulative variables:
% BPV (risky annuity) accumulates the discounted survival probabilities.
BPV = delta_t(1) * survProbs(1) * discounts_CDS(1);
% e accumulates the discounted expected loss (1 - survival probability).
e = discounts_CDS(1) * (1 - survProbs(1));

% Loop over each subsequent CDS maturity to bootstrap survival probabilities.
for i = 2:N
    % Compute the numerator of the survival probability formula:
    % It combines the discounted survival from the previous period with adjustments for 
    % expected loss and risky annuity weighted by the CDS spread.
    num = (1 - recovery) * discounts_CDS(i) * survProbs(i-1) + (1 - recovery) * e - spreadsCDS(i) * BPV;
    
    % Compute the denominator of the formula:
    % This adjusts the discounted survival probability for the spread over the time increment.
    den = ((1 - recovery) + spreadsCDS(i) * delta_t(i)) * discounts_CDS(i);
    
    % Calculate the survival probability for the current period.
    survProbs(i) = num / den;
    
    % Update BPV by adding the current period's contribution: the product of time increment,
    % the new survival probability, and the corresponding discount factor.
    BPV = BPV + delta_t(i) * survProbs(i) * discounts_CDS(i);
    
    % Update the expected loss accumulator 'e' by adding the discounted change in survival probability.
    e = e + discounts_CDS(i) * (survProbs(i-1) - survProbs(i));
    
    % Compute the hazard intensity for the interval [i-1, i]:
    % It is calculated as the logarithmic decrement of survival probability over the time increment.
    intensities(i) = log(survProbs(i-1) / survProbs(i)) / delta_t(i);
end

end
