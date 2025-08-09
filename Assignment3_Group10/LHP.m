function [price] = LHP(pi,Ku,Kd,p,rho)
    % Compute LHP price
    % pi  : recovery of obligor's debt
    % Kd  : amount of losses a portfolio can suffer before tranche's
    %       notional is eroded (TRANCHE INIZIA A PERDERE)
    % Ku  : amount of losses a portfolio can suffer before tranche's is
    %       wiped out (TRANCHE SENZA VALORE)
    % p   : defoult probability of each mortgage
    % rho : correlation betweem mortage

    
    %derivative of inverse of conditional probability
    derivative = @(y) (sqrt(1-rho)/sqrt(rho)) ./ normpdf(norminv(y));
    
    %inverse of conditional probability
    conditionalProbability = @(y) (norminv(y)*sqrt(1-rho)-norminv(p))/sqrt(rho);
    
    %probability density 
    density = @(y) normpdf(conditionalProbability(y))  .* derivative(y);
    
    
    %loss function
    u = Ku/(1-pi);
    d = Kd/(1-pi);
    LossF = @(z) (min(max(z - d,0),(u-d)))/(u-d);
    
    %Integration between 0 and 1
    functionLHP = @(y) LossF(y).*density(y);
    
    LHP = quadgk(functionLHP, 0, 1);
    price = 1 - LHP;
    %LHP2 = sum(functionLHP(z(1:end-1))+functionLHP(z(2:end))+4*functionLHP((z(1:end-1)+z(2:end))/2))*H/6;
    %1 -LHP2
end