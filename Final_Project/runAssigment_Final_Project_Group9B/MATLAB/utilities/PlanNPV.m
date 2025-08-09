function NPV = PlanNPV(AAGR,setDate,flowdates,discounts,c)

zerorates = zeroRates(dates,discounts);
zerorates_flowdates=interp1(dates,zerorates,flowdates,"linear")/100;
B_flowdates = exp(-zerorates_flowdates.*yearfrac(setDate,zerorates_flowdates,3));
C=c*ones(size(flowdates));
for i=12:12:length(C)
    C(i:end)=C(i:end)*AAGR;
end
NPV = C*B_flowdates;


end