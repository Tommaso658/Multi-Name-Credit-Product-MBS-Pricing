function [datesCDS, survProbs, intensities] = bootstrapCDS(datesDF, discounts, datesCDS, spreadsCDS, flag, recovery)
% Computation of the survival probabilities and intensities using
% three different methods based on the 'flag' input.
%
% INPUT:
%   datesDF   : Vector (or array) of dates corresponding to discount factors
%   discounts : Vector (or array) of discount factors
%   datesCDS  : Vector (or array) of CDS maturity/expiry dates
%   spreadsCDS: Vector (or array) of CDS spreads
%   flag      : Choice of the bootstrap method:
%               1 -> Approximation without accrual
%               2 -> With accrual term
%               3 -> Jarrow–Turnbull approach
%   recovery  : Recovery rate for the CDS
%
% OUTPUT:
%   datesCDS   : Same as input, dates of expiry of CDS
%   survProbs  : Computed survival probabilities at each CDS date
%   intensities: Corresponding hazard rates (intensities) for each CDS date

% If the flag is set to 1, call the function that implements
% the bootstrap method without the accrual term.
if flag == 1
    [datesCDS, survProbs, intensities] = bootstrapCDS_NOaccrual( ...
        datesDF, discounts, datesCDS, spreadsCDS, recovery);
end

% If the flag is set to 2, call the function that implements
% the bootstrap method including the accrual term.
if flag == 2
    [datesCDS, survProbs, intensities] = bootstrapCDS_accrual( ...
        datesDF, discounts, datesCDS, spreadsCDS, recovery);
end

% If the flag is set to 3, call the function that implements
% the Jarrow–Turnbull approach (continuously-paid CDS spread).
% Notice here we only use the first discount date (datesDF(1)) 
% because Jarrow–Turnbull might only need an initial discount factor.
if flag == 3
    [datesCDS, survProbs, intensities] = bootstrapCDS_JT( ...
        datesDF(1), datesCDS, spreadsCDS, recovery);
end

end
