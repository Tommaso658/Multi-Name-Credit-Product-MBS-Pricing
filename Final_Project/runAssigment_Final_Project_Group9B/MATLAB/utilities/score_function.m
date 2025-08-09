function score = score_function(HP_prices, copula_prices_up, copula_prices_down)
% SCORE_FUNCTION computes the proportion of HP_prices that lie between 
% the lower and upper bounds for each tranche (column-wise).
%
% Inputs:
%   - HP_prices: Matrix (n x m) of observed prices, with each column corresponding to a tranche.
%   - copula_prices_up: Matrix (n x m) or scalar/vector of upper bounds.
%   - copula_prices_down: Matrix (n x m) or scalar/vector of lower bounds.
%
% Output:
%   - score: Row vector (1 x m), each value is the proportion of HP_prices(:,i) within the bounds.

    % Ensure input matrices have consistent size
    [n, ~] = size(HP_prices);

    % Logical matrix of values within bounds
    within_bounds = (HP_prices > copula_prices_down) & (HP_prices < copula_prices_up);

    % Compute the proportion of values within bounds for each tranche (column)
    score = sum(within_bounds, 1) ./ n; % 1 x m row vector

end