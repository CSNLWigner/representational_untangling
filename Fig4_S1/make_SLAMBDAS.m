% pre-generating and saving wavelength (spatial period) samples
% from its lognormal distribution

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

M1 = 100000; M2 = 200000;

pmu = 0.95;
psigma = 0.55;
lognorm = @(x) exp(-(log(x)-pmu).^2/psigma^2/2)./x/psigma/sqrt(2*pi);

lambda_range = [0.5 7.5];

s1seed = 1;
s2seed = 2;

tic

SLAMBDAS1 = samples(lambda_range(1),lambda_range(2),lognorm,M1,'permutated_representative',s1seed);
SLAMBDAS2 = samples(lambda_range(1),lambda_range(2),lognorm,M2,'permutated_representative',s2seed);

runtime = toc

save('SLAMBDAS_data.mat');