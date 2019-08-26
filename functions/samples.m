function X = samples(xmin,xmax,density_function,N,method,random_seed)
% used functions: rejection_sampler, bdd_splitter
% X: column vector of generated samples
% xmin: lower bound of the domain of the density function
% xmax: upper bound of the domain of the density function
% density_function: function handle of the unnormalized positive 1D density
% N: number of samples (N >> 1)
% method: sampling method describing how samples are generated
%  'independent': independent samples
%  'representative': representative deterministic samples
%  'bin-uniform': uniform within non-overlapping representative bins
% random_seed: random seed value for non-deterministic methods (optional)
switch method
    case 'independent'
        if exist('random_seed')
            X = rejection_sampler(xmin,xmax,density_function,N,random_seed);
        else
            X = rejection_sampler(xmin,xmax,density_function,N);
        end
    case 'representative'
        X1 = bdd_splitter(xmin,xmax,density_function,N);
        Y = arrayfun(density_function,X1);
        X = ((X1(1:N).*Y(1:N)+X1(2:N+1).*Y(2:N+1))./(Y(1:N)+Y(2:N+1)))';
        if exist('random_seed')
            rng(random_seed);
            X = X(randperm(N));
        end
    case 'permutated_representative'
        X1 = bdd_splitter(xmin,xmax,density_function,N);
        Y = arrayfun(density_function,X1);
        X = ((X1(1:N).*Y(1:N)+X1(2:N+1).*Y(2:N+1))./(Y(1:N)+Y(2:N+1)))';
        if exist('random_seed')
            rng(random_seed);
        end
        X = X(randperm(N));
    case 'bin-uniform'
        if exist('random_seed')
            rng(random_seed);
        end
        X1 = bdd_splitter(xmin,xmax,density_function,N)';
        noise = rand(N,1);
        noise(noise==1) = 0; % guarantee no overlapping
        X = X1(1:N)+noise.*diff(X1);
    otherwise
        error('Unexpected method!');
end