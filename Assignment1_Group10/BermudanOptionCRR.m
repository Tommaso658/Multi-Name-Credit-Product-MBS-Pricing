function optionPrice = BermudanOptionCRR(F0,K,B,T,sigma,N,flag)
%BermudanOptionCRR calcola il prezzo di un'opzione bermudiana 
% utilizzando il metodo ad albero CRR.
%
% L'opzione ha maturità T (in anni) ed è esercitabile anticipatamente
% ogni 1 mese. Nel caso specifico:
%   T = 4/12 (4 mesi) e N = 80 (step dell'albero),
%   dunque ogni esercizio intermedio avviene ogni (1/12)/(T/N) = 20 step.
%
% INPUT
%   F0   : forward price (prezzo iniziale)
%   K    : strike price
%   B    : discount factor su T, B = exp(-r*T)
%   div  : dividend yield (non utilizzato in questo esempio)
%   T    : time-to-maturity (in anni) --> per 4 mesi, T = 4/12
%   sigma: volatilità
%   N    : numero di step (qui N = 80)
%   flag : 1 per call, -1 per put
%
% OUTPUT
%   optionPrice: prezzo dell'opzione bermudiana

% Calcolo del passo temporale
delta_t = T/N;

% Calcolo dei parametri CRR
delta_x = sigma * sqrt(delta_t);
u = exp(delta_x);
d = exp(-delta_x);    % d = 1/u
q = (1 - d) / (u - d);
r = -log(B)/T;        % da B = exp(-r*T)
Bfrac = exp(-r * delta_t);  % fattore di sconto per ogni passo

% Calcolo dell'intervallo in step per l'opportunità di early exercise:
% 1 mese = 1/12 di anno, quindi:
exercise_interval = round((1/12) / delta_t);  % in questo caso 20

% Costruzione dell'albero dei prezzi del sottostante.
% Tree(i+1,j+1) rappresenta il prezzo al tempo j*delta_t con i down moves.
for j = 0:N
    for i = 0:j
        Tree(i+1,j+1) = F0 * u^(j - 2*i);
    end
end

% Inizializzazione dell'albero dei prezzi dell'opzione.
PriceTree = Tree;

% Calcolo del payoff terminale (al tempo T).
if flag == 1   % Call
    PriceTree(:,end) = max(Tree(:,end) - K, 0);
else           % Put
    PriceTree(:,end) = max(K - Tree(:,end), 0);
end

% Backward induction per il calcolo del prezzo dell'opzione.
for j = N:-1:1
    for i = 1:j
        PriceTree(i,j) = Bfrac * (q * PriceTree(i,j+1) + (1 - q) * PriceTree(i+1,j+1));
    end
    
    % Se il tempo corrisponde a un'opportunità di early exercise (ogni 1 mese)
    % allora aggiorniamo il valore al nodo confrontando con il payoff immediato.
    if mod(j, exercise_interval) == 0
        if flag == 1   % Call: esercizio = S - K
            for i = 1:j
                PriceTree(i,j) = max(Tree(i,j) - K, PriceTree(i,j));
            end
        else           % Put: esercizio = K - S
            for i = 1:j
                PriceTree(i,j) = max(K - Tree(i,j), PriceTree(i,j));
            end
        end
    end
end

% Il prezzo dell'opzione è il valore alla radice dell'albero.
optionPrice = PriceTree(1,1);
end
