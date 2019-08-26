% script for generating data for Fig. 3
% generated data: performance curve of the linear decoder: no nuisance
% extrenal functions used:
%   paramset
%   randisk
%   cresponse
%   nonlin
%   perfeval_fc
%   perfeval_pfc
%   maxfit

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

% see Table 1 in the paper (MK: M_theta, disk_radius: R, V0: u_DC, V1: u_AC)
N = 500; K = 10; MK = 10; M1 = 2000; M2 = 5000;

maxiter = 15000; conv_crit = 1e-3; eta0 = 1e-9;

disk_radius = 90; gcenter_seed = 1; gphi_seed = 1;

stim_DC = 0.5; V0 = -60; V1 = 12; noise = 3; 

lambda0 = 3; sigma0 = 2; c0 = 0.5;

phi0 = 180;

thresholds = -75:0.5:-40;

N1 = length(thresholds);

FC_LIN = zeros(1,N1);
PFC_LIN = zeros(1,N1);
RUNTIMES_LIN = zeros(1,N1);
ITERATIONS = zeros(1,N1);
W_DATA = zeros(N1,N+1,K-1);

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
s0phi = {'const',phi0};
s0DC = {'const',stim_DC};
s0AC = {'const',stim_DC*c0};

SBANK0 = paramset(s0lambda,s0theta,s0phi,s0DC,s0AC);
SBANK1 = repmat(SBANK0,M1,1);
SBANK2 = repmat(SBANK0,M2,1);

% ORIENTATION CATEGORIES:

Y1 = ceil(SBANK1(:,2)/(180/K));
Y2 = ceil(SBANK2(:,2)/(180/K));

% MP RESPONSE:

U_MEAN = V0+V1*cresponse(GBANK,SBANK0);

U1 = repmat(U_MEAN,M1,1)+noise*randn(K*MK*M1,N);
U2 = repmat(U_MEAN,M2,1)+noise*randn(K*MK*M2,N);

W0 = 0;
for i1 = 1:N1

    threshold = thresholds(i1);
    
    tic
    
    R1 = nonlin(U1,threshold,1);
    R2 = nonlin(U2,threshold,1);
    
    [W1,etas] = mnrfitbb(R1,Y1,W0,maxiter,conv_crit,eta0);
    iterations = length(etas);
    if iterations == 0
        [W1,etas] = mnrfitbb(R1,Y1,0,maxiter,conv_crit,eta0);
        iterations = length(etas); 
    end
    ITERATIONS(i1) = iterations;
    ETAS(i1).vals = etas;
    W_DATA(i1,:,:) = W1;
    W0 = W1;
    
    if i1 == 1 fprintf('%d',iterations); else fprintf('-%d',iterations); end
    
    PK_LIN = mnrvals(W1,R2);
    
    FC_LIN(i1) = perfeval_fc(PK_LIN,Y2);
    PFC_LIN(i1) = perfeval_pfc(PK_LIN,Y2);
    
    RUNTIMES_LIN(i1) = toc;
end

runtime = ceil(sum(RUNTIMES_LIN(:))/60);
minute_data = mod(runtime,60);
hour_data = (runtime-minute_data)/60;
fprintf(' %dh%dm\n',hour_data,minute_data);

clear SBANK0 SBANK1 SBANK2 Y1 Y2 U_MEAN U1 U2 R1 R2 W0 W1 PK_LIN

save(sprintf('%s_data.mat',mfilename));

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
plot([])
hold on
axis([thresholds(1) thresholds(end) 0 1])
xlabel('threshold [mV]');
ylabel('performance');
title({sprintf('N = %d, K = %d, MK = %d, M1 = %d, M2 = %d',N,K,MK,M1,M2),sprintf('convcrit = %g, eta0 = %g, resolution = %g [mV] (runtime: %d)',conv_crit,eta0,thresholds(2)-thresholds(1),runtime)})
p1 = plot(thresholds,FC_LIN,'color','k','linewidth',2);
p2 = plot(thresholds,PFC_LIN,'color','k','linewidth',2,'linestyle',':');
l0 = line([thresholds(1) thresholds(end)],[1/K 1/K],'color','r','linestyle',':');
[maxx_fc,maxy_fc] = maxfit(thresholds,FC_LIN,3);
[maxx_pfc,maxy_pfc] = maxfit(thresholds,PFC_LIN,3);
line([maxx_fc maxx_fc],[0 0.033],'color','k','linewidth',2)
line([maxx_pfc maxx_pfc],[0 0.033],'color','k','linewidth',2,'linestyle',':')
scatter(maxx_fc,maxy_fc,100,'r','d','filled')
scatter(maxx_pfc,maxy_pfc,100,'r','d','filled')
legend([p1,p2,l0],{'FC','PFC','chance level'})

saveas(gcf,sprintf('%s_done.png',mfilename));