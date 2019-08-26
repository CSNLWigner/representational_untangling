function PK = mnrvals(W,X)
% vectorized quick mnrval applied to several x vectors simultaneously
% W is an output weight matrix comes from mnrfit: (N+1) x (K-1)
% X contains feature vectors: M x N
% PK contains class probabilities for each feature vector: K x M
[M,N] = size(X);
K = size(W,2)+1;
B0 = [W';zeros(1,N+1)];
X1 = [ones(1,M);X'];
E = exp(B0*X1);
S = sum(E);
PK = E./repmat(S,K,1);
end