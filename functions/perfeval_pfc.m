function p = perfeval_pfc(PK,Y)
% gives back probabilistic fraction correct
% PK: posterior probabilities for categories
% Y: real categories (from 1 to K)
M = length(Y);
S = 0;
for m = 1:M
    S = S+log(PK(Y(m),m));
end
p = exp(S/M);
end