function [integr] = probabilityMExact(m, I, p, rho)
    %Compute the probability density of having m defoults in HP
    % m   : number of obligors defoulted in the reference portfolio
    % I   : total obligors in the reference portfolio
    % p   : defoult probability of each mortgage
    % rho : correlation betweem mortage

    %computation of the binomial coefficient
    warning('off', 'all')
    binCoef = nchoosek(I,m);
    
    % defining the functions 
    % conditional probability
    conditionalProbability = @(y) normcdf((norminv(p)-sqrt(rho)*y)/sqrt(1-rho));
    % integrand of ptobability density
    integrand = @(y) normpdf(y)*binCoef.*(conditionalProbability(y).^m).*((1-conditionalProbability(y)).^(I-m));
    
    % integral for obtain probability density
    integr  = quadgk(integrand, -6,6);
end
