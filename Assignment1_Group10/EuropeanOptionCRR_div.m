function optionPrice = EuropeanOptionCRR_div(F0,K,B,div,T,sigma,N,flag)
% INPUT:
%   F0   : prezzo iniziale del sottostante (spot)
%   K    : strike
%   B    : discount factor su T, B = exp(-r*T)
%   div  : dividend yield continuo (es. 0.05 -> 5%)
%   T    : time-to-maturity (in anni)
%   sigma: volatilità
%   N    : numero di step binomiali
%   flag : 1 per call, -1 per put
%
% OUTPUT:
%   optionPrice : prezzo dell'opzione EU w div

% 1) Passo temporale
delta_t = T / N;

% 2) Calcolo del tasso risk-free r (da B = exp(-r*T))
r = -log(B) / T;

% 3) Parametri CRR (up, down) e probabilità risk-neutral con dividendo
u = exp(sigma * sqrt(delta_t));
d = exp(-sigma * sqrt(delta_t));
% <-- Qui la modifica per includere 'div' -->
q = (exp((r - div)*delta_t) - d) / (u - d);

% 4) Fattore di sconto per ciascun passo
discount = exp(-r * delta_t);



% 6) Albero dei prezzi del sottostante
Tree = zeros(N+1, N+1);
for j = 0:N
    for i = 0:j
        % j-i up-moves, i down-moves
        Tree(i+1, j+1) = F0 * (u^(j - i)) * (d^i);
    end
end

% 7) Albero dei payoff dell'opzione
PriceTree = zeros(N+1, N+1);

% Payoff terminale a scadenza (j=N)
if flag == 1  % Call
    PriceTree(:, end) = max(Tree(:, end) - K, 0);
else          % Put
    PriceTree(:, end) = max(K - Tree(:, end), 0);
end

% 8) Backward induction
for j = N:-1:1
   
    for i = 1:j
        PriceTree(i, j) = discount * ...
            ( q * PriceTree(i, j+1) + (1 - q) * PriceTree(i+1, j+1) );
    end
    
end

% 9) Prezzo dell'opzione = valore alla radice dell'albero
optionPrice = PriceTree(1, 1);

end
