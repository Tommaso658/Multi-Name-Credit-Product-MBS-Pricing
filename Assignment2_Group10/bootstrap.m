function [dates, discounts]=bootstrap(datesSet, ratesSet)


%% building a complete set of swaps dates

%creating a vector of suitable length
oldSwapsLength = size(datesSet.swaps,1);
newSwapsLength = round((datesSet.swaps(oldSwapsLength,1)-datesSet.swaps(1,1))/365.25)+1;
newSwapsDates = zeros(newSwapsLength,1); 

%for loop to compute a complete set of swaps' expiry dates
index = 1; 
for year = 2025:2025+newSwapsLength-1 %this 
   newSwapsDates(index,1) = datenum(year,2,1);
   index = index +1;
end 

%going to the following business day in case of holiday
newSwapsDates = busdate(newSwapsDates); 


%% building a complete set of swaps rates 
newSwapsRates = zeros(newSwapsLength,1);
indexNew = 1; %inizializing index for loop
indexOld = 1;

%this loops inserts the rates found in the market in the new swapsRates table
%and sets at NaN the missing rates
for index = 1:newSwapsLength 
    
    if datesSet.swaps(indexOld,1) == newSwapsDates(indexNew,1)
        newSwapsRates(indexNew,1) = (ratesSet.swaps(indexOld,1)+ratesSet.swaps(indexOld,2))/2;
        indexOld = indexOld+1;
    else 
        newSwapsRates(indexNew,1) = NaN;
    end 
    indexNew = indexNew+1;

end 

%spline interpolation for the rates we don't have
missingIndices = isnan(newSwapsRates);
newSwapsRates(missingIndices,1) = interp1(yearfrac(datesSet.settlement, newSwapsDates(~missingIndices,1),3), newSwapsRates(~missingIndices,1), yearfrac(datesSet.settlement, newSwapsDates(missingIndices,1),3), 'spline');

%updating the struct and creating the vectors for discount factors and
%output dates
datesSet.swaps = newSwapsDates;
ratesSet.swaps = newSwapsRates;
discounts = zeros(newSwapsLength+4+4,1);
dates = zeros(newSwapsLength+4+4,1);


%% Bootstrap with deposits

%using known formulas to get discount factors B from quoted rates
for j = 1:4 
    dates(j) = datesSet.depos(j,1);
    Lrate = (ratesSet.depos(j,1)+ratesSet.depos(j,2))/2;
    discounts(j) = 1/(1+(datesSet.depos(j)-datesSet.settlement)*Lrate/360);
end 


%% Bootstrap with futures 
gap = 4; %gap to equate index inside futures table with index in general date and discount table

%this for loop distinguishes the cases when we need extrapolation
%and when we need interpolation
%the computes the rates given by the quoted futures and the 
%desired bootstrap rate
for j = 1:7
   Lrate = (ratesSet.futures(j,1)+ ratesSet.futures(j,2))/2;

   if dates(j+gap-1)-datesSet.futures(j,1)>0 %if the value date of the future is before the last known discount's date
                                          %we interpolate using linear zero
                                          %rate

       %we compute the zero rates then interpolate                                   
       yKnownt1 = -log(discounts(j+gap-1))*365/(dates(j+gap-1)-datesSet.settlement);
       yKnownt0 = -log(discounts(j+gap-2))*365/(dates(j+gap-2)-datesSet.settlement);
       yInterp = interp1(dates(j+gap-2:j+gap-1),[yKnownt0, yKnownt1], datesSet.futures(j,1), 'linear');
       discountInterp = exp(-(datesSet.futures(j,1)-datesSet.settlement)*yInterp/365);
       discountAtValueDate = discountInterp;
   end 

   if dates(j+gap-1) <= datesSet.futures(j)%here we do flat extrapolation or nothing 
       discountAtValueDate = discounts(gap+j-1);
   end 

   time = datesSet.futures(j,2)-datesSet.futures(j,1);
   time = time/360;
   forwardDiscount = 1/(1+time*Lrate); %compute the forward discount
   %multiply the forward discount for the discount at value date to obtain the discount at expiry
   discounts(j+gap)=discountAtValueDate*forwardDiscount; 
   dates(j+gap) = datesSet.futures(j,2);
end 

%% Bootstrap with swaps
gap = 8; %adjusting the gap to connect the internal swap index to the general one 

%interpolation to get the discount for the first swap coupon

%computing zero rates
yKnownt1 = -log(discounts(gap))*365/(dates(gap)-datesSet.settlement);
yKnownt0 = -log(discounts(gap-1))*365/(dates(gap-1)-datesSet.settlement);
%interpolating
yInterp = interp1(dates(gap-1:gap),[yKnownt0, yKnownt1], datesSet.swaps(1,1)-365.25, 'linear');
%going back to discount
swapDiscountOne = exp(-(datesSet.swaps(1,1)-365.25-datesSet.settlement)*yInterp/365);

%we define new vectors to make the loop easier
swapsDatesForLoop = [datesSet.settlement, datesSet.swaps(1,1)-365.25,  datesSet.swaps'];
swapsDiscount = [swapDiscountOne, zeros(newSwapsLength,1)'];
gap = 11;%adjusting the gap to connect the internal swap index to the general one 
%this loop does the bootstrap 
for j = 1:newSwapsLength
   discountSum = 0;
   for i = 1:j
       discountSum = discountSum + swapsDiscount(i)*yearfrac(swapsDatesForLoop(i),swapsDatesForLoop(i+1),6);
   end 
   swapsDiscount(j+1) = (1-ratesSet.swaps(j,1)*discountSum)/(1+ratesSet.swaps(j,1)*yearfrac(swapsDatesForLoop(j+1),swapsDatesForLoop(j+2),6));
   dates(gap+j) = datesSet.swaps(j);
end 

%% final vectors, plotting and tables
discounts = [1 discounts(1:gap)' swapsDiscount(2:end)];
dates = [datesSet.settlement dates'];
end



