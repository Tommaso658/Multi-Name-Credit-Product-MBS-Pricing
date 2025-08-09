function settlement_date=change_convention_date(date)

%the goal is to change the date of settlement if 2 days after it is in the weekend

if isbusday(date)==1
    settlement_date=date;
else
    settlement_date = date;
    while isbusday(settlement_date)==0
        settlement_date = settlement_date-1;
    end
end
