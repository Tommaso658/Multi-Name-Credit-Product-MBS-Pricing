function [prob,intensity] = bootstrap_CDS(datesDF, discounts, datesCDS, spreadsCDS, flag, recovery)
    %The goal is to bootstrap the intensity supposed piecewise constant
    settlement = datesCDS(1);
    new_discounts = getDiscount(discounts,datesCDS,datesDF,settlement);
    act365 = 3;
    prob = 1;
    prob = [prob,(1-recovery)/((1-recovery)+spreadsCDS(1)*yearfrac(settlement,datesCDS(2),act365))];
    intensity = -log(prob(end))/yearfrac(settlement,datesCDS(2),act365);
    for i=2:(length(datesCDS)-1)
        num_new_prob_1 = (1-recovery)*(sum(discounts(1:i-1).*(prob(1:i-1)-prob(2:i)))+discounts(i)*prob(i));
        num_new_prob_2 = spreadsCDS(i)*sum(yearfrac(datesCDS(1:i-1),datesCDS(2:i),act365).*prob(2:i).*discounts(1:i-1));
        num_new_prob = num_new_prob_1-num_new_prob_2;
        den_new_prob = yearfrac(datesCDS(i),datesCDS(i+1),act365)*discounts(i)+(1-recovery)*discounts(i);
        new_prob=num_new_prob./den_new_prob;
        prob = [prob, new_prob];
        new_lamb = -(log(prob(end)/prob(end-1)))/yearfrac(datesCDS(i),datesCDS(i+1),act365);
        intensity = [intensity,new_lamb];
    end
end
