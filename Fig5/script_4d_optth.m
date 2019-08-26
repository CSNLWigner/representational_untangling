% script for generating data (optimal thresholds and decoder weights) for script_4d_lognorm.m

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

N = 500; K = 10; MK = 10; J = 50;

M1 = 200000; M2 = 500000;

maxiter = 10000; conv_crit = 1e-3; eta0 = 1e-8;

disk_radius = 90; gcenter_seed = 1; gphi_seed = 1;

glambda = 3; gsigma = 2;

stim_DC = 0.5; V0 = -60; V1 = 12; noise = 3; 

thresholds = -58.4:0.2:-55.2;

N1 = length(thresholds);

pmu = 0.95;
psigma = 0.55;
lognorm = @(x) exp(-(log(x)-pmu).^2/psigma^2/2)./x/psigma/sqrt(2*pi);

palpha = 2.4;
pbeta = 3.6;
beta = @(x) (x.^(palpha-1)).*(1-x).^(pbeta-1);

nuisance(1).name = 'phase';
nuisance(1).range = [0 360];

nuisance(2).name = 'wavelength';
nuisance(2).range = [0.5 7.5];

nuisance(3).name = 'contrast';
nuisance(3).range = [0 1];

FC_4D = zeros(1,N1);
PFC_4D = zeros(1,N1);
ITERATIONS_4D = zeros(1,N1);
RUNTIMES_4D = zeros(1,N1);

% GABOR POPULATION:

gset_x0 = {'const',0};
gset_y0 = {'const',0};
gset_lambda = {'const',glambda};
gset_theta = {'grid',0,180,N};
gset_phi = {'perm',0,360,gphi_seed};
gset_sigma = {'const',gsigma};
GBANK = paramset(gset_x0,gset_y0,gset_lambda,gset_theta,gset_phi,gset_sigma);
GBANK(:,1:2) = randisk(disk_radius,N,gcenter_seed);

s0lambda = {'const',0};
s0theta = {'grid',0,180,K*MK};
s0phi = {'grid',0,360,J};
s0DC = {'const',stim_DC};
s0AC = {'const',stim_DC*0};

tic

s1seed = 1;
s2seed = 2;
SBANK0 = paramset(s0lambda,s0theta,s0phi,s0DC,s0AC);
SBANK1 = repmat(SBANK0,M1/size(SBANK0,1),1);
SBANK2 = repmat(SBANK0,M2/size(SBANK0,1),1);
SBANK1(:,1) = samples(nuisance(2).range(1),nuisance(2).range(2),lognorm,M1,'permutated_representative',s1seed);
SBANK1(:,5) = stim_DC*samples(0,1,beta,M1,'permutated_representative',s1seed+5);
SBANK2(:,1) = samples(nuisance(2).range(1),nuisance(2).range(2),lognorm,M2,'permutated_representative',s2seed);
SBANK2(:,5) = stim_DC*samples(0,1,beta,M2,'permutated_representative',s2seed+5);

save(sprintf('%s_data.mat',mfilename));

Y1 = ceil(SBANK1(:,2)/(180/K));
Y2 = ceil(SBANK2(:,2)/(180/K));

U1 = V0+V1*cresponse(GBANK,SBANK1)+noise*randn(M1,N);
U2 = V0+V1*cresponse(GBANK,SBANK2)+noise*randn(M2,N);

toc

W0 = 0;
for i1 = 1:N1
    
    threshold = thresholds(i1);
    
    tic;
    
    R1 = nonlin(U1,threshold,1);
    R2 = nonlin(U2,threshold,1);
    
    [W1,etas] = mnrfitbb(R1,Y1,W0,maxiter,conv_crit,eta0);
    iterations = length(etas);
    W0 = W1;
    
    ITERATIONS_4D(i1) = iterations;
    if i1 == 1 fprintf('%d',iterations); else fprintf('-%d',iterations); end
    
    PK = mnrvals(W1,R2);
    
    FC_4D(i1) = perfeval_fc(PK,Y2);
    PFC_4D(i1) = perfeval_pfc(PK,Y2);
    
    RUNTIMES_4D(i1) = toc/60; % minutes
end
fprintf('\n');

subplot(2,1,1)

hold on
plot(thresholds,FC_4D,'color',[.8 .8 .8],'linewidth',2);
axis([thresholds(1) thresholds(end) 0.9*min(FC_4D) 1.1*max(FC_4D)])
xlabel('threshold [mV]');
ylabel('fraction correct');
title(sprintf('N = %d, K = %d, MK = %d, M1 = %d, M2 = %d',N,K,MK,M1,M2))

[maxx,maxy,xfit,yfit] = maxfit(thresholds,FC_4D,5,100);
scatter(maxx,maxy,20,'k','d','filled')
plot(xfit,yfit,'color','r','linewidth',1);
line([maxx maxx],[0.9*min(FC_4D) 1.1*max(FC_4D)],'linestyle',':')

OPTTH_4D = maxx;

subplot(2,1,2)

hold on
plot(thresholds,PFC_4D,'color',[.8 .8 .8],'linewidth',2);
axis([thresholds(1) thresholds(end) 0.9*max(PFC_4D) 1.1*max(PFC_4D)])
xlabel('threshold [mV]');
ylabel('probabilistic fraction correct');

[maxx,maxy,xfit,yfit] = maxfit(thresholds,PFC_4D,5,100);
scatter(maxx,maxy,20,'k','d','filled')
plot(xfit,yfit,'color','r','linewidth',1);
line([maxx maxx],[0.9*min(PFC_4D) 1.1*max(PFC_4D)],'linestyle',':')

saveas(gcf,sprintf('%s_done.png',mfilename));

R1 = nonlin(U1,OPTTH_4D,1);
[W4D,etas] = mnrfitbb(R1,Y1,0,maxiter,conv_crit,eta0);
 
clear SBANK0 SBANK1 SBANK2 Y1 Y2 U1 U2 R1 R2 W0 W1 PK
save(sprintf('%s_data.mat',mfilename));