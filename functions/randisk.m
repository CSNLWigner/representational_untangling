function XY = randisk(R,N,random_seed)
% sampling from 2D uniform distribution on a disk
% XY: Nx2 output matrix with random coordinate pairs
% R: disk radius (default = 1)
% N: number of independent samples (default = 1)
% random_seed: random seed value (optional)
if nargin < 1
    R = 1;
end
if nargin < 2
    N = 1;
end
if exist('random_seed')
    rng(random_seed);
end
XY = zeros(N,2);
for n = 1:N
    while 1 % rejection sampling
        x = 2*R*rand()-R;
        y = 2*R*rand()-R;
        if x^2+y^2 <= R^2, break, end
    end
    XY(n,:) = [x y];
end
end