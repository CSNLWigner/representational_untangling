% script for generating data for Fig. 5

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

% load optimal thresholds and weights of the optimal linear decoder
% calculated with script_4d_optth.m
load('script_4d_optth_data.mat','OPTTH_4D','W4D')

N = 500; K = 10; MK = 10; J = 50;

M1 = 200000; M2 = 500000;

maxiter = 10000; conv_crit = 1e-3; eta0 = 1e-8;

disk_radius = 90; gcenter_seed = 1; gphi_seed = 1;

glambda = 3; gsigma = 2;

stim_DC = 0.5; V0 = -60; V1 = 12; noise = 3; 

nuisance(1).name = 'phase';
nuisance(1).label = '$$\varphi$$';
nuisance(1).range = [0 360];

nuisance(2).name = 'wavelength';
nuisance(2).label = '$$\lambda$$';
nuisance(2).range = [0.01 20];

nuisance(3).name = 'contrast';
nuisance(3).label = '$$c$$';
nuisance(3).range = [0 1];

palpha = 2.4;
pbeta = 3.6;
beta = @(x) (x.^(palpha-1)).*(1-x).^(pbeta-1);

pmu_values = [0.95 0.5 1.4 1.85];
psigma_values = [0.55 0.55 0.55 0.55];

thresholds = -58.8:0.2:-54.8;

N1 = length(pmu_values);
N2 = length(thresholds);

FC_4D = zeros(N1,N2);
ITERATIONS_4D = zeros(N1,N2);
RUNTIMES_4D = zeros(N1,N2);

OPT_TH = zeros(N1,1);
OPT_FC = zeros(N1,1);
GREY_FC = zeros(N1,1);

% GABOR POPULATION:

gset_x0 = {'const',0};
gset_y0 = {'const',0};
gset_lambda = {'const',glambda};
gset_theta = {'grid',0,180,N};
gset_phi = {'perm',0,360,gphi_seed};
gset_sigma = {'const',gsigma};
GBANK = paramset(gset_x0,gset_y0,gset_lambda,gset_theta,gset_phi,gset_sigma);
GBANK(:,1:2) = randisk(disk_radius,N,gcenter_seed);

tic

s0lambda = {'const',0};
s0theta = {'grid',0,180,K*MK};
s0phi = {'grid',0,360,J};
s0DC = {'const',stim_DC};
s0AC = {'const',stim_DC*0};
SBANK0 = paramset(s0lambda,s0theta,s0phi,s0DC,s0AC);
SBANK1 = repmat(SBANK0,M1/size(SBANK0,1),1);
SBANK2 = repmat(SBANK0,M2/size(SBANK0,1),1);
s1seed = 1;
s2seed = 2;
SBANK1(:,5) = stim_DC*samples(0,1,beta,M1,'permutated_representative',s1seed+5);
SBANK2(:,5) = stim_DC*samples(0,1,beta,M2,'permutated_representative',s2seed+5);

toc

for i1 = 1:N1
    
    pmu = pmu_values(i1);
    psigma = psigma_values(i1);
    lognorm = @(x) exp(-(log(x)-pmu).^2/psigma^2/2)./x/psigma/sqrt(2*pi);
    
    SBANK1(:,1) = samples(nuisance(2).range(1),nuisance(2).range(2),lognorm,M1,'permutated_representative',s1seed);
    SBANK2(:,1) = samples(nuisance(2).range(1),nuisance(2).range(2),lognorm,M2,'permutated_representative',s2seed);
    
    Y1 = ceil(SBANK1(:,2)/(180/K));
    Y2 = ceil(SBANK2(:,2)/(180/K));
    
    U1 = V0+V1*cresponse(GBANK,SBANK1)+noise*randn(M1,N);
    U2 = V0+V1*cresponse(GBANK,SBANK2)+noise*randn(M2,N);
    
W0 = 0;
for i2 = 1:N2
    
    tic;
    
    threshold = thresholds(i2);
    
    R1 = nonlin(U1,threshold,1);
    R2 = nonlin(U2,threshold,1);
    
    [W1,etas] = mnrfitbb(R1,Y1,W0,maxiter,conv_crit,eta0);
    iterations = length(etas);
    W0 = W1;
    
    ITERATIONS_4D(i1,i2) = iterations;
    if i2 == 1 fprintf('%d',iterations); else fprintf('-%d',iterations); end
    
    PK = mnrvals(W1,R2);
    
    FC_4D(i1,i2) = perfeval_fc(PK,Y2);
    
    RUNTIMES_4D(i1,i2) = toc/60; % minutes
end

[maxx,maxy] = maxfit(thresholds,FC_4D(i1,:),5);
OPT_TH(i1) = maxx;
OPT_FC(i1) = maxy;

R2 = nonlin(U2,OPTTH_4D,1);
PK = mnrvals(W4D,R2);
GREY_FC(i1) = perfeval_fc(PK,Y2);

fprintf('\n');

end

clear SBANK0 SBANK1 SBANK2 Y1 Y2 U1 U2 R1 R2 W0 W1 PK
save(sprintf('%s_data.mat',mfilename));

plot([])
hold on
xlabel('FRNL threshold [mV]')
ylabel('fraction correct')
xlim([thresholds(1) thresholds(end)])
ylim([0.95*(1/K) 1.1*max(FC_4D(:))])
line([thresholds(1) thresholds(end)],[1/K 1/K],'linestyle',':','color','k')    

PLOTS = 1:N1;
LEGENDS = {};
for i1 = 1:N1
    color = [0 1-i1/N1/2 0];
    PLOTS(i1) = plot(thresholds,FC_4D(i1,:),'linewidth',2,'color',color);
    line([OPT_TH(i1) OPT_TH(i1)],[1/K OPT_FC(i1)],'linewidth',2,'color',color)
    line([thresholds(1) thresholds(end)],[GREY_FC(i1) GREY_FC(i1)],'linewidth',1,'color',color)    
    LEGENDS{end+1} = sprintf('mu = %g, sig = %g',pmu_values(i1),psigma_values(i1));
end
legend(PLOTS,LEGENDS)
legend boxoff

saveas(gcf,sprintf('%s_done.png',mfilename));