function loss = lossFunction(i,pi,Ku,Kd,p,rho);
    % Compute the Loss fuction
    % i   : total obligors in the reference portfolio
    % pi  : recovery of obligor's debt
    % Kd  : amount of losses a portfolio can suffer before tranche's
    %       notional is eroded (TRANCHE INIZIA A PERDERE)
    % Ku  : amount of losses a portfolio can suffer before tranche's is
    %       wiped out (TRANCHE SENZA VALORE)
    % p   : defoult probability of each mortgage
    % rho : correlation betweem mortage

    u = Ku/(1-pi);
    d = Kd/(1-pi);
      
    totalLoss = 0;
    
    %calculation of probability density and loss fuction for every possibility of obligors defoulted
    for m = 0:i
        integral = probabilityMExact(m, i, p, rho);   %probability density
        z = m/i;
        loss = (min(max(z - d,0),(u-d)))/(u-d);   %loss function
        totalLoss = loss*integral+totalLoss;   %I sum loss for otain total loss
    end
    
    %setup the output
    loss = totalLoss;

end