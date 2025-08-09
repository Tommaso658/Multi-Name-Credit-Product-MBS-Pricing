function optionPrice = BermudanOptionCRR_div(F0,K,B,div,T,sigma,N,flag)
%BermudanOptionCRR calcola il prezzo di un'opzione bermudiana 
% (call o put) con albero binomiale CRR, tenendo conto di un
% dividend yield continuo 'div' nel calcolo della probabilità risk-neutral.
%
% L'opzione è esercitabile alla fine di ogni mese, fino a maturità T.
%
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
%   optionPrice : prezzo dell'opzione bermudiana

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

% 5) Calcolo dell'intervallo in step per esercizio mensile
%    (esempio: se T=4/12 e N=80, exercise_interval=20)
exercise_interval = round((1/12) / delta_t);

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
    % Valore di continuazione
    for i = 1:j
        PriceTree(i, j) = discount * ...
            ( q * PriceTree(i, j+1) + (1 - q) * PriceTree(i+1, j+1) );
    end
    
    % Se siamo a una data mensile, confrontiamo con payoff di esercizio
    if mod(j, exercise_interval) == 0 && j > 0 && j < N
        for i = 1:j
            if flag == 1  % Call
                immediatePayoff = max(Tree(i, j) - K, 0);
            else          % Put
                immediatePayoff = max(K - Tree(i, j), 0);
            end
            PriceTree(i, j) = max(immediatePayoff, PriceTree(i, j));
        end
    end
end

% 9) Prezzo dell'opzione = valore alla radice dell'albero
optionPrice = PriceTree(1, 1);

end
