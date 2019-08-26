% pre-generating and saving contrast samples
% from its beta distribution

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

stim_DC = 0.5;

M1 = 100000; M2 = 200000;

palpha = 2.4;
pbeta = 3.6;
beta = @(x) (x.^(palpha-1)).*(1-x).^(pbeta-1);

s1seed = 1;
s2seed = 2;

tic

SCONTRASTS1 = stim_DC*samples(0,1,beta,M1,'permutated_representative',s1seed+5);
SCONTRASTS2 = stim_DC*samples(0,1,beta,M2,'permutated_representative',s2seed+5);

runtime = toc

save('SCONTRASTS_data.mat');