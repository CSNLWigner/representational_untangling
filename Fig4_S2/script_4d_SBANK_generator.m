% stimulus bank pre-generator

close all
clearvars

N = 500; K = 10; MK = 10; J = 50;

M1 = 100000; M2 = 200000;

nuisance(1).name = 'phase';
nuisance(1).label = '$$\varphi$$';
nuisance(1).range = [0 360];
nuisance(1).mean = 180;

nuisance(2).name = 'wavelength';
nuisance(2).label = '$$\lambda$$';
nuisance(2).range = [0.5 7.5];
nuisance(2).mean = 3;

nuisance(3).name = 'contrast';
nuisance(3).label = '$$c$$';
nuisance(3).range = [0 1];
nuisance(3).mean = 0.4;

phi0 = 180;
lambda0 = 3;
c0 = 0.5; 

pmu = 0.95;
psigma = 0.55;
lognorm = @(x) exp(-(log(x)-pmu).^2/psigma^2/2)./x/psigma/sqrt(2*pi);

palpha = 2.4;
pbeta = 3.6;
beta = @(x) (x.^(palpha-1)).*(1-x).^(pbeta-1);

stim_DC = 0.5;

s0lambda = {'const',lambda0};
s0theta = {'grid',0,180,K*MK};
s0phi = {'grid',0,360,J};
s0DC = {'const',stim_DC};
s0AC = {'const',stim_DC*c0};

SBANK0 = paramset(s0lambda,s0theta,s0phi,s0DC,s0AC);
SBANK1 = repmat(SBANK0,M1/size(SBANK0,1),1);
SBANK2 = repmat(SBANK0,M2/size(SBANK0,1),1);

tic

s1seed = 1;
s2seed = 2;
SBANK1(:,1) = samples(nuisance(2).range(1),nuisance(2).range(2),lognorm,M1,'permutated_representative',s1seed);
SBANK1(:,5) = stim_DC*samples(0,1,beta,M1,'permutated_representative',s1seed+5);
SBANK2(:,1) = samples(nuisance(2).range(1),nuisance(2).range(2),lognorm,M2,'permutated_representative',s2seed);
SBANK2(:,5) = stim_DC*samples(0,1,beta,M2,'permutated_representative',s2seed+5);

toc

cd(fileparts(which(mfilename)));
save('SBANK.mat','SBANK1','SBANK2');