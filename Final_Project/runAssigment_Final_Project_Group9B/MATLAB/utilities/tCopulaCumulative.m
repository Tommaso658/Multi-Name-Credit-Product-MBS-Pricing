function [priceCum, priceCumUp, priceCumDown, varExpLoss] = tCopulaCumulative(discount, Kd, KuList, Nsim, rho, recovery, pd, I, nu)
% --------------------------------------------------------------------------------------------------------
% tCopulaCumulative
%
% Computes cumulative tranche prices using the t-Student copula model via Monte Carlo simulation.
% This model accounts for fat tails and tail dependence, unlike the Gaussian copula.
%
% INPUTS:
%   discount   – discount factor for present value calculation
%   Kd         – lower attachment point (shared by all tranches)
%   KuList     – vector of upper attachment points (each defines a tranche)
%   Nsim       – number of Monte Carlo simulations
%   rho        – asset correlation (copula parameter)
%   recovery   – recovery rate on defaults
%   pd         – marginal default probability (over the time horizon)
%   I          – number of exposures in the portfolio
%   nu         – degrees of freedom of the t-distribution
%
% OUTPUTS:
%   priceCum      – vector of expected tranche prices (% of notional)
%   priceCumUp    – upper bound of 1σ confidence interval (price + std error)
%   priceCumDown  – lower bound of 1σ confidence interval (price - std error)
%   varExpLoss    – variance of the loss estimator (for each tranche)
% --------------------------------------------------------------------------------------------------------

% Step 1 – Compute the default threshold K from the inverse t-distribution
K  = tinv(pd, nu);   % Scalar threshold such that P(X < K) = pd

%% 1. Generate t-distributed latent factors (systemic and idiosyncratic)

% Systemic risk factor y ~ t_nu (using standard normal scaled by chi-squared)
y  = randn(Nsim,1) ./ sqrt(chi2rnd(nu, Nsim, 1) / nu);  % Nsim × 1
% size(y)  % Just for debug/inspection

% Idiosyncratic components (unique for each of the I names)
eps = randn(Nsim, I) ./ sqrt(chi2rnd(nu, Nsim, I) / nu);  % Nsim × I
% size(eps)  % Just for debug/inspection

% Latent variables: combine systemic and idiosyncratic risk
X = sqrt(rho) .* y + sqrt(1 - rho) .* eps;  % Nsim × I
% Each row is a scenario of the I latent variables

%% 2. Compute defaults and total portfolio loss for each scenario

% Default occurs when X_ij < K (threshold)
Default = X < K;  % Nsim × I logical array (true if default)

% Total portfolio loss (sum of defaulted exposures × LGD)
L = (1 - recovery) / I .* sum(Default, 2);  % Nsim × 1 loss vector

%% 3. Compute cumulative tranche loss between Kd and Ku (per scenario)

% Ensure KuList is a row vector
Ku = KuList(:)';  % 1 × m (m = number of tranches)

% Numerator: loss falling into each tranche [Kd, Ku]
numer = max(L - Kd, 0) - max(L - Ku, 0);  % Nsim × m

% Denominator: width of each tranche
denom = Ku - Kd;  % 1 × m

% Tranche payoff as a percentage of the tranche notional
payoff = numer ./ denom;  % Nsim × m

%% 4. Estimate expected tranche loss and variance

% Expected loss (per tranche)
expLoss = mean(payoff, 1);  % 1 × m

% Variance of the Monte Carlo estimator (not the payoff variance itself)
varExpLoss = var(payoff, 0, 1) / Nsim;  % 1 × m

% Standard error of the expected loss
stdExpLoss = 1.96*sqrt(varExpLoss);  % 1 × m

%% 5. Compute discounted tranche prices and 1σ confidence intervals

% Base case: price = 1 - expected loss
priceCum = discount .* (1 - expLoss);  % 1 × m

% Upper and lower bounds with ±1 standard error
priceCumUp   = discount .* (1 - (expLoss - stdExpLoss));  % optimistic estimate
priceCumDown = discount .* (1 - (expLoss + stdExpLoss));  % conservative estimate
end
