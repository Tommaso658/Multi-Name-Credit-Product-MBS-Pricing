function [priceCum, priceCumUp, priceCumDown, varExpLoss] = gaussian_copula(discount, Kd, KuList, Nsim, rho, recovery, pd, I)
% -------------------------------------------------------------------------
% Gaussian one-factor copula – Monte-Carlo pricing of cumulative tranches
%
% Inputs
%   discount  : discount factor (PV of 1 € at maturity)
%   Kd        : lower attachment of the equity tranche (usually 0)
%   KuList    : row-vector of upper detachments   e.g. [0.03 0.06 0.09]
%   Nsim      : number of Monte-Carlo scenarios   (e.g. 1e5)
%   rho       : asset correlation ρ  (scalar in (0,1))
%   recovery  : recovery rate R  (LGD = 1−R)
%   pd        : unconditional default probability Q over the horizon
%   I         : portfolio size (number of identical loans)
%
% Outputs
%   priceCum     : discounted cumulative price for each 0–Ku(i) (1×m)
%   priceCumUp   : upper 1-σ Monte-Carlo bound on priceCum        (1×m)
%   priceCumDown : lower 1-σ Monte-Carlo bound on priceCum        (1×m)
%   varExpLoss   : variance of the Monte-Carlo mean of the loss   (1×m)
% -------------------------------------------------------------------------

%% 1.  Generate latent factors
K      = norminv(pd);            % default threshold in normal space
Y      = randn(Nsim,1);          % systemic factor      (N×1)
Eps    = randn(Nsim,I);          % idiosyncratic shocks (N×I)
X      = sqrt(rho).*Y + sqrt(1-rho).*Eps;

%% 2.  Default indicators and portfolio loss
Default = X < K;                                        % N×I logical
L       = (1-recovery)/I .* sum(Default,2);             % N×1 loss

%% 3.  Cumulative tranche pay-offs 0–Ku
Ku      = KuList(:)';                                   % ensure 1×m
numer   = max(L-Kd,0) - max(L-Ku,0);                    % N×m (implicit expansion)
denom   = Ku - Kd;                                      % 1×m
payoff  = numer ./ denom;                               % N×m relative loss

%% 4.  Monte-Carlo mean and variance
expLoss     = mean(payoff,1);                           % 1×m
varExpLoss  = var(payoff,0,1) / Nsim;                   % variance of the estimator
stdExpLoss  = 1.96*sqrt(varExpLoss);                         % 1×m standard error

%% 5.  Discounted cumulative prices and error bounds
priceCum      = discount .* (1 - expLoss);
priceCumUp    = discount .* (1 - (expLoss - stdExpLoss));  % +1 σ
priceCumDown  = discount .* (1 - (expLoss + stdExpLoss));  % −1 σ
end