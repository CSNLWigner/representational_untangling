% script for generating data for Fig. 4 figure supplement 1
% data: decoding performance from membrane potentials (MP)
% extrenal functions used:
%   paramset
%   randisk
%   cresponse
%   nonlin
%   perfeval_fc
%   perfeval_pfc
% used pre-generated data: 'SLAMBDAS_data.mat', 'SCONTRASTS_data.mat'

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

N = 500; K = 10; MK = 10; J = 50;

M1 = 100000; M2 = 200000;

maxiter = 10000; conv_crit = 1e-3; eta0 = 1e-8;

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

% defining lognormal distribution for wavelength (spatial period)
pmu = 0.95;
psigma = 0.55;
lognorm = @(x) exp(-(log(x)-pmu).^2/psigma^2/2)./x/psigma/sqrt(2*pi);

% defining beta distribution for contrast
palpha = 2.4;
pbeta = 3.6;
beta = @(x) (x.^(palpha-1)).*(1-x).^(pbeta-1);

disk_radius = 90;

gseed_values = 1:25;

sigma0 = 2; stim_DC = 0.5; V0 = -60; V1 = 12; noise = 3; 

threshold = -80;

N1 = length(gseed_values);
N2 = 1+3+3+1; 

FC_MP = zeros(N1,N2);
PFC_MP = zeros(N1,N2);
ITERATIONS = zeros(N1,N2);
RUNTIMES = zeros(N1,N2);

load('SLAMBDAS_data.mat','SLAMBDAS1','SLAMBDAS2')
load('SCONTRASTS_data.mat','SCONTRASTS1','SCONTRASTS2')
     
for i1 = 1:N1
    
gphi_seed = gseed_values(i1);
gcenter_seed = 100+gseed_values(i1);

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

for i2 = 1:N2
    
    tic;
    
    s0lambda = {'const',lambda0};
    s0theta = {'grid',0,180,K*MK};
    s0phi = {'const',phi0};
    s0DC = {'const',stim_DC};
    s0AC = {'const',stim_DC*c0};
    
    if (i2 == 2) || (i2 == 5) || (i2 == 6) || (i2 == 8)
        s0phi = {'grid',0,360,J};
    end
    
    SBANK0 = paramset(s0lambda,s0theta,s0phi,s0DC,s0AC);
    SBANK1 = repmat(SBANK0,M1/size(SBANK0,1),1);
    SBANK2 = repmat(SBANK0,M2/size(SBANK0,1),1);
    
    s1seed = 1;
    s2seed = 2;
    switch i2
        case {3, 5}
            SBANK1(:,1) = SLAMBDAS1;
            SBANK2(:,1) = SLAMBDAS2;
        case {4, 6}
            SBANK1(:,5) = SCONTRASTS1;
            SBANK2(:,5) = SCONTRASTS2;
        case {7, 8}
            SBANK1(:,1) = SLAMBDAS1;
            SBANK1(:,5) = SCONTRASTS1;
            SBANK2(:,1) = SLAMBDAS2;
            SBANK2(:,5) = SCONTRASTS2;
    end
    
    Y1 = ceil(SBANK1(:,2)/(180/K));
    Y2 = ceil(SBANK2(:,2)/(180/K));
    
    U1 = V0+V1*cresponse(GBANK,SBANK1)+noise*randn(M1,N);
    U2 = V0+V1*cresponse(GBANK,SBANK2)+noise*randn(M2,N);
    
    R1 = nonlin(U1,threshold,1);
    R2 = nonlin(U2,threshold,1);
    
    [W1,etas] = mnrfitbb(R1,Y1,0,maxiter,conv_crit,eta0);
    ITERATIONS(i1,i2) = length(etas);
    
    PK = mnrvals(W1,R2);
    
    FC_MP(i1,i2) = perfeval_fc(PK,Y2);
    PFC_MP(i1,i2) = perfeval_pfc(PK,Y2);
    
    RUNTIMES(i2,i1) = toc/60; % minutes
    
    fprintf('.');
end

fprintf('\n');

end

clear SBANK0 SBANK1 SBANK2 Y1 Y2 U1 U2 R1 R2 PK
save([mfilename,'_data.mat']);