function X = rejection_sampler(xmin,xmax,density_function,N,random_seed)
% X: column vector of generated samples
% xmin: lower bound of the domain of the density function
% xmax: upper bound of the domain of the density function
% density_function: function handle of the unnormalized positive 1D density
% N: number of independent samples
% random_seed: random seed value (optional)
domain_width = xmax-xmin;
X = NaN*zeros(N,1);
if exist('random_seed')
    rng(random_seed);
end
x0 = fminbnd(@(x) -density_function(x),xmin,xmax);
upper_bound = 1.1*density_function(x0);
for n = 1:N
    while 1 % rejection sampling
        x_value = xmin+domain_width*rand();
        y_value = upper_bound*rand();
        f_value = density_function(x_value);
        if f_value < 0
            error('Density must be non-negative!');
        end
        if y_value <= f_value, break, end
    end
    X(n) = x_value;
end
end