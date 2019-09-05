% script for generating data for Fig. 8
% extrenal functions used:
%   randsphere 
%   paramset
%   randisk
%   cresponse
%   nonlin
%   mnrfitbb
%   mnrvals
%   perfeval_fc

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

load('heterogeneous_matrixdata.mat')

[KAPPAMATRIX,NOISEMATRIX] = meshgrid(kappa_mvalues,noise_mvalues);

cell_noise = repmat(NOISEMATRIX(:),20,1)';
cell_kappa = repmat(KAPPAMATRIX(:),20,1)';

OPTTH2 = OPTTH';
cell_optth = repmat(OPTTH2(:),20,1)';

rmax = 7;
N_points = 100;
dim = 500;
R_values = rmax/N_points:rmax/N_points:rmax;
X = randsphere(N_points,dim,1);
POINTS = [];
for n = 1:N_points
    point = X(n,:);
    radius = sqrt(sum(point.^2));
    POINTS = [POINTS; R_values(n)*point/radius];
end
POINTS = [zeros(1,dim);POINTS];
R_values = [0,R_values];

N = 500; K = 10; MK = 10; J = 50; M1 = 25; M2 = 50;

maxiter = 5000; conv_crit = 1e-4; eta0 = 1e-8;

disk_radius = 90; gcenter_seed = 1; gphi_seed = 1;

stim_DC = 0.5; V0 = -60; V1 = 12; prefactor = 16.7;

lambda0 = 3; sigma0 = 2; c0 = 0.5;

PERFORMANCES = zeros(1,N_points+1);
RUNTIMES = zeros(1,N_points+1);
ITERATIONS = zeros(1,N_points+1);

% GABOR POPULATION:
  
gx0 = {'const',0};
gy0 = {'const',0};
glambda = {'const',lambda0};
gtheta = {'grid',0,180,N};
gphi = {'perm',0,360,gphi_seed};
gsigma = {'const',sigma0};
GBANK = paramset(gx0,gy0,glambda,gtheta,gphi,gsigma);
GBANK(:,1:2) = randisk(disk_radius,N,gcenter_seed);

% STIMULUS BANK:

s0lambda = {'const',lambda0};
s0theta = {'grid',0,180,K*MK};
s0phi = {'grid',0,360,J};
s0DC = {'const',stim_DC};
s0AC = {'const',stim_DC*c0};

SBANK0 = paramset(s0lambda,s0theta,s0phi,s0DC,s0AC);
SBANK1 = repmat(SBANK0,M1,1);
SBANK2 = repmat(SBANK0,M2,1);

Y1 = ceil(SBANK1(:,2)/(180/K));
Y2 = ceil(SBANK2(:,2)/(180/K));

U_MEAN = V0+V1*cresponse(GBANK,SBANK0);

NOISE1 = repmat(cell_noise,K*MK*J*M1,1);
NOISE2 = repmat(cell_noise,K*MK*J*M2,1);

U1 = repmat(U_MEAN,M1,1)+NOISE1.*randn(K*MK*J*M1,N);
U2 = repmat(U_MEAN,M2,1)+NOISE2.*randn(K*MK*J*M2,N);

W0 = 0;
for n = 1:N_points+1
    
    point = POINTS(n,:);
    cell_thresholds = V0+(1-point).*(cell_optth-V0);

    R1 = U1;
    R2 = U2;
    for nn = 1:N
        R1(:,nn) = nonlin(U1(:,nn),cell_thresholds(nn),prefactor,cell_kappa(nn));
        R2(:,nn) = nonlin(U2(:,nn),cell_thresholds(nn),prefactor,cell_kappa(nn));
    end
    
    tic
    
    [W1,etas] = mnrfitbb(R1,Y1,W0,maxiter,conv_crit,eta0);
    ITERATIONS(n) = length(etas);
    W0 = W1;
    
    PK = mnrvals(W1,R2);
    
    PERFORMANCES(n) = perfeval_fc(PK,Y2);
    
    RUNTIMES(n) = toc;
    
    runtime = ceil(RUNTIMES(n)/60);
    minute_data = mod(runtime,60);
    hour_data = (runtime-minute_data)/60;

    fprintf('iterations: %d, runtime: %dh%dm\n',ITERATIONS(n),hour_data,minute_data)
end

clear SBANK0 SBANK1 SBANK2 Y1 Y2 U_MEAN U1 U2 R1 R2 W0 W1 PK

save(sprintf('%s_data.mat',mfilename));