function p = perfeval_fc(PK,Y)
% gives back fraction correct with MAP
% PK: posterior probabilities for categories
% Y: real categories (from 1 to K)
[~,max_indices] = max(PK);
D = max_indices(:)-Y(:);
p = sum(~D)/length(Y);
end