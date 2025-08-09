function [DV01, BPV, DV01_z] = sensSwap(setDate,fixedLegPaymentDates,fixedRate,dates,discounts,discounts_DV01)
%the function computes sensitivities: DV01, BPV, DV01_z
%
%INPUTS:

%setDate:                           Settlement Date
%fixedLegPaymentDates:              Payment Dates of the Fixed Leg
%fixedRate:                         Swap Rate (fixed)
%dates:                             Dates of the discount factors bootstrap
%discounts:                         Bootstrap discount factors
%discounts_DV01:                    Bootstrap discount with 1bp shift

%OUTPUTS:

%DV01                               NPV_shifted(t0) - NPV(t0)
%BPV                                Basis Point Value
%DV01_z                             Diff NPV when shifted rates 


% Deltas between payments using Act/365 convention ->option 6
delta_t = yearfrac([setDate;fixedLegPaymentDates(1:end-1)],fixedLegPaymentDates,6);

% Convert discounts (obtained from bootstrap) into zeroRates
z_Rates_boot = zeroRates(dates,discounts); 
z_Rates_boot = z_Rates_boot/100;

z_Rates_DV01= zeroRates(dates,discounts_DV01); 
z_Rates_DV01 = z_Rates_DV01/100;

% data convertion = Act/365 
yf_fixedleg = yearfrac(setDate, fixedLegPaymentDates, 3);
yf_dates = yearfrac(setDate, dates, 3);


% interpolation to compute the missing zRates
z_Rates_fix = interp1(yf_dates, z_Rates_boot, yf_fixedleg,'linear','extrap');
z_Rates_fixDV01 = interp1(yf_dates,z_Rates_DV01,yf_fixedleg,'linear','extrap');

% move from rates to discounts
B = exp(-z_Rates_fix.*yf_fixedleg);
B_DV01 = exp(-z_Rates_fixDV01.*yf_fixedleg);

% computation of the net present value (NPV)
BPV = sum(delta_t.*B);
BPV = BPV*(1e-4);
NPV_fixed = sum(delta_t.*B)*fixedRate;
NPV_float = 1 - B(end);
NPV = NPV_fixed - NPV_float;

%BPV_shifted = sum(delta_t.*B_DV01);
NPV_fixed_shifted = sum(delta_t.*B_DV01)*fixedRate;
NPV_float_shifted = 1 - B_DV01(end);
NPV_shifted = NPV_fixed_shifted - NPV_float_shifted;

% DV01(t0) = NPV_shifted(t0) - NPV(t0)
DV01 = abs(NPV_shifted-NPV);


% Shift of the zeroRates of 1bps
zRates_boot_shifted = z_Rates_boot + 1e-4;

% Interpolation to find the desired zeroRates on the shifted curve
zRates_Shifted = interp1(dates(2:end),zRates_boot_shifted(2:end),fixedLegPaymentDates,'linear','extrap'); 
% move from rates to discounts
B_z = exp(-zRates_Shifted.*yf_fixedleg);

%computation of the net present value(NPV)
BPV_z = sum(delta_t.*B_z);
NPV_fixed_z = fixedRate*BPV_z;
NPV_float_z = 1 - B_z(end);
NPV_z = NPV_fixed_z - NPV_float_z;

% DV01_z(t0) = NPV_z(t0) - NPV(t0)
DV01_z = abs(NPV - NPV_z);

end