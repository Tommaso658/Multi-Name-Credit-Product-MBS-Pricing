function optionPrice = EuropeanOptionKOCRR(F0,K,KO,B,T,sigma,N)
%European option UP&OUT price with CRR tree approach
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility


delta_t = T/N;

%parameters of the CRR model
delta_x = sigma*sqrt(delta_t);
u = exp(delta_x);
d = exp(-delta_x); %d = 1/u
q = (1 - d)/(u-d);
r = -1/T * log(B);
Bfrac = exp(-r * delta_t);

%we create the empty tree
Tree = zeros(N+1,N+1);

for j=0:N
    for i=0:j
        Tree(i+1,j+1)=F0*u^(-2*i+j);
    end
end

%we are pricing a call option
Tree(:,end) = max(Tree(:,end)-K,0).*(Tree(:,end)<KO);


%going backward in time
for j=N:-1:1
    for i = 1:j
        Tree(i,j) = Bfrac*(q*Tree(i,j+1) + (1-q)*Tree(i+1,j+1));
    end
end

%actualizing the value at the root of the tree
optionPrice = Tree(1,1);
end

