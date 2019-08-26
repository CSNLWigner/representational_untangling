% script for generating data for Fig. 4
% generated data: performance curves for all nuisance combinations
% extrenal functions used:
%   paramset
%   randisk
%   samples
%   cresponse
%   nonlin
%   mnrfitbb
%   mnrvals
%   perfeval_fc
%   perfeval_pfc

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

% see Table 1 in the paper (MK: M_theta, disk_radius: R, V0: u_DC, V1: u_AC)
N = 500; K = 10; MK = 10; J = 50;

M1 = 100000; M2 = 200000; % size of the stimulus bank for training (M1) and for test (M2)

maxiter = 10000; conv_crit = 1e-3; eta0 = 1e-8;

% possible nuisance parameters and their settings: range, etc.
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

% distribution for wavelength (spatial period)
pmu = 0.95;
psigma = 0.55;
lognorm = @(x) exp(-(log(x)-pmu).^2/psigma^2/2)./x/psigma/sqrt(2*pi);

% distribution for contrast
palpha = 2.4;
pbeta = 3.6;
beta = @(x) (x.^(palpha-1)).*(1-x).^(pbeta-1);

disk_radius = 90; gcenter_seed = 1; gphi_seed = 1;

sigma0 = 2; stim_DC = 0.5; V0 = -60; V1 = 12; noise = 3; 

thresholds = -75:0.5:-40;

N1 = 1+3+3+1; N2 = length(thresholds);

FC = zeros(N1,N2);
PFC = zeros(N1,N2);
ITERATIONS = zeros(N1,N2);
RUNTIMES = zeros(N1,N2);

% GABOR POPULATION:

gx0 = {'const',0};
gy0 = {'const',0};
glambda = {'const',lambda0};
gtheta = {'grid',0,180,N};
gphi = {'perm',0,360,gphi_seed};
gsigma = {'const',sigma0};
GBANK = paramset(gx0,gy0,glambda,gtheta,gphi,gsigma);
GBANK(:,1:2) = randisk(disk_radius,N,gcenter_seed);

warning off

for i1 = 1:N1
    
    fprintf('nuisance combination %d/%d:\n',i1,N1);
    
    s0lambda = {'const',lambda0};
    s0theta = {'grid',0,180,K*MK};
    s0phi = {'const',phi0};
    s0DC = {'const',stim_DC};
    s0AC = {'const',stim_DC*c0};
    
    if (i1 == 2) || (i1 == 5) || (i1 == 6) || (i1 == 8)
        s0phi = {'grid',0,360,J};
    end
    
    SBANK0 = paramset(s0lambda,s0theta,s0phi,s0DC,s0AC);
    SBANK1 = repmat(SBANK0,M1/size(SBANK0,1),1);
    SBANK2 = repmat(SBANK0,M2/size(SBANK0,1),1);
    
    s1seed = 1;
    s2seed = 2;
    switch i1
        case {3, 5}
            SBANK1(:,1) = samples(nuisance(2).range(1),nuisance(2).range(2),lognorm,M1,'permutated_representative',s1seed);
            SBANK2(:,1) = samples(nuisance(2).range(1),nuisance(2).range(2),lognorm,M2,'permutated_representative',s2seed);
        case {4, 6}
            SBANK1(:,5) = stim_DC*samples(0,1,beta,M1,'permutated_representative',s1seed+5);
            SBANK2(:,5) = stim_DC*samples(0,1,beta,M2,'permutated_representative',s2seed+5);
        case {7, 8}
            SBANK1(:,1) = samples(nuisance(2).range(1),nuisance(2).range(2),lognorm,M1,'permutated_representative',s1seed);
            SBANK1(:,5) = stim_DC*samples(0,1,beta,M1,'permutated_representative',s1seed+5);
            SBANK2(:,1) = samples(nuisance(2).range(1),nuisance(2).range(2),lognorm,M2,'permutated_representative',s2seed);
            SBANK2(:,5) = stim_DC*samples(0,1,beta,M2,'permutated_representative',s2seed+5);
    end
    
    Y1 = ceil(SBANK1(:,2)/(180/K));
    Y2 = ceil(SBANK2(:,2)/(180/K));
    
    U1 = V0+V1*cresponse(GBANK,SBANK1)+noise*randn(M1,N);
    U2 = V0+V1*cresponse(GBANK,SBANK2)+noise*randn(M2,N);
    
W0 = 0;
for i2 = 1:N2
    
    threshold = thresholds(i2);
    
    tic;
    
    R1 = nonlin(U1,threshold,1);
    R2 = nonlin(U2,threshold,1);
    
    [W1,etas] = mnrfitbb(R1,Y1,W0,maxiter,conv_crit,eta0);
    iterations = length(etas);
    ITERATIONS(i1,i2) = iterations;
    W0 = W1;
    
    if i2 == 1 fprintf('%d',iterations); else fprintf('-%d',iterations); end
    
    PK = mnrvals(W1,R2);
    
    FC(i1,i2) = perfeval_fc(PK,Y2);
    PFC(i1,i2) = perfeval_pfc(PK,Y2);
    
    RUNTIMES(i1,i2) = toc/60;
end

runtime_sum = ceil(sum(RUNTIMES(i1,:)));
minute_data = mod(runtime_sum,60);
hour_data = (runtime_sum-minute_data)/60;
fprintf(' %dh%dm\n',hour_data,minute_data);

end

clear SBANK0 SBANK1 SBANK2 Y1 Y2 U1 U2 R1 R2 PK
save([mfilename,'_data.mat']);