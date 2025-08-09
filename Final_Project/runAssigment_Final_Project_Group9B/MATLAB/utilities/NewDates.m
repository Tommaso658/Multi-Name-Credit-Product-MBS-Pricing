function [new_dates,new_swap_rates] = NewDates(dates,swap_rates)
    new_dates=(dates(1)-1):365:dates(end);
    for i=1:length(new_dates)
        if leapyear(new_dates(i))==1
            new_dates(i)=new_dates(i)+1;
        end
    end
    new_swap_rates = zeros(size(new_dates));
    for i=2:length(new_dates)
        new_dates(i)=change_convention_date(new_dates(i));
    end
    new_swap_rates = spline(dates,swap_rates,new_dates);
end