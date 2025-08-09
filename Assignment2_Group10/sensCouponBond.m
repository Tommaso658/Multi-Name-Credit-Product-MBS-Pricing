function MacD = sensCouponBond(setDate, couponPaymentDates, fixedRate, dates, discounts)
% The function computes the Macaulay Duration of a coupon bond 

%INPUTS:

% setDate:              Settlement Date
% couponPaymentDates:   Payments Dates of the CB
% fixedRate:            ZCB rate
% dates:                Dates of the DF computed in bootstrap
% discounts:            Discount Factors from bootstrap

%OUTPUT:

% MacD:                 Macaulay Duration

%Compute deltas
deltas=yearfrac(setDate.*ones(length(couponPaymentDates),1),couponPaymentDates,6);

%Compute discounts at Payment Dates
zRates = zeroRates(dates, discounts);
B=interp1(dates(2:end),zRates(2:end),couponPaymentDates,'linear','extrap');
B=exp(-B/100.*deltas);

%Compute price of Coupon Bond
coupons=0.*couponPaymentDates;
coupons=fixedRate+coupons;coupons(end)=coupons(end)+1;
Price=sum(coupons.*B);

%Compute Duration
MacD=sum(coupons.*deltas.*B)/Price;

end


