function X = bdd_splitter(xmin,xmax,density_function,N,reltol)
% bdd_splitter: bounded density domain splitter
%  goal: subdivision of the domain into bins of equal mass
% X: row vector of boundary points of bins in order
%  X(1) = xmin, X(end) = xmax, length(X) = number of bins + 1
% xmin: lower bound of the domain of the density function
% xmax: upper bound of the domain of the density function
% density_function: function handle of the (unnormalized) density
%  density: non-negative 1D function having a monotone comulative function
% N: number of bins
% reltol (optional): relative resolution for numeric integrals
if nargin < 5
    reltol = 0.001;
end
total_mass = integral(density_function,xmin,xmax,'reltol',reltol/N);
bin_mass = total_mass/N;
% adaptive bisection to reach the necessary resolution:
XX = [xmin,xmax];
AA = [total_mass];
while 1
    KK = find(AA > 0.9*bin_mass);
    for i = 1:length(KK)
        k = KK(i);
        XX = [XX(1:k),(XX(k)+XX(k+1))/2,XX(k+1:end)];
        AA = [AA(1:k),AA(k:end)];
        AA(k) = integral(density_function,XX(k),XX(k+1),'reltol',reltol);
        AA(k+1) = integral(density_function,XX(k+1),XX(k+2),'reltol',reltol);
        KK = KK+1;
    end
    if isempty(KK)
        break;
    end
end
% end of adaptive bisection subrutine
SS = [0,cumsum(AA)];
SS = (total_mass/SS(end))*SS; % correction
X = [xmin,zeros(1,N-1),xmax];
k = 0;
for n = 1:N-1
    cumulative_mass = n*bin_mass;
    k = k+find(SS(k+1:end) > cumulative_mass,1);
    f1 = density_function(XX(k-1));
    f2 = density_function(XX(k));
    slope = (f2-f1)/(XX(k)-XX(k-1));    
    X(n+1) = XX(k-1)+(sqrt(f1^2+2*slope*(cumulative_mass-SS(k-1)))-f1)/slope;
end
end