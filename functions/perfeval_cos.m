function p = perfeval_cos(PK,Y)
% gives back cosine error
% PK: posterior probabilities for categories
% Y: real categories (from 1 to K)
M = length(Y);
K = max(Y);
S = 0;
for m = 1:M
    [~,maxk] = max(PK(:,m));
    dtheta = (maxk-Y(m))*pi/K; 
    S = S+cos(2*dtheta);
end
p = S/M;
end