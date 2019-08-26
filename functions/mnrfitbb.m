function [W,etas] = mnrfitbb(X,Y,W0,maxiter,convcrit,eta0,logdisp)
% Barzalai-Borwein method to perform gradient descent (Barzilai and Borwein, 1988)
% also see: https://en.wikipedia.org/wiki/Gradient_descent
% input X: M x N
% input Y: M x 1 /values: 1,2,...,K/
% input W0: (1+N) x (K-1)
% output W: (1+N) x (K-1)
N = size(X,2);
K = max(Y);
M = length(Y);
if nargin < 7, logdisp = false; end
if nargin < 6, eta0 = 1e-5; end
if nargin < 5, convcrit = 1e-3; end
if nargin < 4, maxiter = 100; end
if nargin < 3, W0 = zeros(1+N,K-1); end
if W0 == 0, W0 = zeros(1+N,K-1); end
X1 = [ones(M,1),X]; % M x (1+N)
T = zeros(M,K);
for k = 1:K
    T(Y==k,k) = 1;
end
T1 = T(:,1:K-1); % M x (K-1)
eta1 = eta0;
W = W0;
etas = [];
E12 = zeros(1,2);
iter = 0;
while iter < maxiter
    iter = iter + 1;
    P1 = exp(X1*W); % M x (K-1)
    sumcol = sum(P1,2)+1; % M x 1
    P1 = P1./repmat(sumcol,1,K-1);
    if iter <= 2
        P = [P1 1-sum(P1,2)];
        E12(iter) = -sum(log(P(T==1))); % energy
    end
    if (iter == 2) && (E12(2) >= E12(1))
        eta1 = eta1/5;
        W = W0;
        etas = [];
        if eta1 < 1e-15, break; end
        E12 = zeros(1,2);
        iter = 0;
        continue;
    end
    if iter ~= 1, DE_old = DE; end
    DE = X1'*(P1-T1); % (1+N) x (K-1)
    if iter == 1
        eta = eta1;
    else
        DDE = DE-DE_old;
        DWDDE = sum(sum((W-W_old).*(DDE)));
        DDE_norm2 = sum(sum(DDE.^2));
        eta = DWDDE/DDE_norm2;
        if (isinf(eta)) || (eta == 0 ) || (isnan(eta))
            if logdisp, fprintf('etabreak\n'); end
            break;
        end
    end
    etas(end+1) = eta;
    W_old = W;
    W = W-eta*DE;
    if (iter > 1) && (~any(abs(W(:)-W_old(:)) > convcrit*abs(W_old(:)))), break; end
end
end