function [datesCDS, survProbs, intensities] = bootstrapCDS_JT(SettlementDate, datesCDS, spreadsCDS, recovery)
% bootstrapCDS_JT applies the Jarrow–Turnbull (JT) method to calculate survival 
% probabilities and hazard rates (intensities) for a set of CDS contracts.
%
% INPUT:
%   SettlementDate : The current date (valuation date)
%   datesCDS       : Vector of CDS expiry (maturity) dates
%   spreadsCDS     : Vector of CDS spreads
%   recovery       : Recovery rate assumed for the CDS
%
% OUTPUT:
%   datesCDS   : The same input vector of CDS expiry dates
%   survProbs  : Calculated survival probabilities at each CDS expiry date
%   intensities: Hazard rates (intensities) corresponding to each CDS expiry date

% -------------------------------------------------------------------------
% Calculate the hazard rates (intensities) using the JT method.
% In this method, the intensity is assumed constant over time and is given by:
%       intensity = spread / (1 - recovery)
%
% Here, spreadsCDS is divided by (1 - recovery) to obtain the hazard rate.
intensities = spreadsCDS / (1 - recovery);

% -------------------------------------------------------------------------
% Compute the survival probabilities assuming the hazard rate (lambda) is 
% constant over each year. First, convert the CDS expiry dates into year 
% fractions relative to the SettlementDate using day-count convention 3.
dates_CDS_yf = yearfrac(SettlementDate, datesCDS, 3);
delta_t = diff([0;dates_CDS_yf]);

% Then, compute the survival probability at each CDS expiry date using the 
% exponential decay formula:
%       survProb = exp( - intensity * time )
% where 'time' is the year fraction from the SettlementDate.
survProbs = exp(-intensities .* dates_CDS_yf);

% intensities assuming lambda as averages of the step functions lambda
intensities = -log(survProbs./[1;survProbs(1:end-1)])./delta_t; 

end
