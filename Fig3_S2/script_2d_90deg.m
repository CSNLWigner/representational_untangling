% script for creating data for Fig. 3 figure supplement 2
% created data: 25 performance curves (R = 90deg, phase nuisance)
% extrenal functions used:
%   paramset
%   randisk
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
N = 500; K = 10; MK = 10; J = 50; M1 = 10; M2 = 20;

maxiter = 5000; conv_crit = 1e-4; eta0 = 1e-8;

disk_radius = 90; gcenter_seed = 1; gphi_seed = 1;

stim_DC = 0.5; V0 = -60; V1 = 12; noise = 3; 

lambda0 = 3; sigma0 = 2; c0 = 0.5;

thresholds = -75:0.5:-40;

N1 = 25;
N2 = length(thresholds);

FC_LIN = zeros(N1,N2);
PFC_LIN = zeros(N1,N2);
RUNTIMES = zeros(N1,N2);
ITERATIONS = zeros(N1,N2);

for i1 = 1:N1
    
    % GABOR POPULATION:

    gx0 = {'const',0};
    gy0 = {'const',0};
    glambda = {'const',lambda0};
    gtheta = {'grid',0,180,N};
    gphi = {'perm',0,360,gphi_seed};
    gsigma = {'const',sigma0};
    GBANK = paramset(gx0,gy0,glambda,gtheta,gphi,gsigma);
    GBANK(:,1:2) = randisk(disk_radius,N,i1);

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

    rng(i1)
    U1 = repmat(U_MEAN,M1,1)+noise*randn(K*MK*J*M1,N);
    U2 = repmat(U_MEAN,M2,1)+noise*randn(K*MK*J*M2,N);

    W0 = 0;
    for i2 = 1:N2

        threshold = thresholds(i2);

        tic

        R1 = nonlin(U1,threshold,1);
        R2 = nonlin(U2,threshold,1);

        [W1,etas] = mnrfitbb(R1,Y1,W0,maxiter,conv_crit,eta0);
        ITERATIONS(i1,i2) = length(etas);
        W0 = W1;

        if i2 == 1 fprintf('%d',ITERATIONS(i1,i2)); else fprintf('-%d',ITERATIONS(i1,i2)); end

        PK = mnrvals(W1,R2);

        FC_LIN(i1,i2) = perfeval_fc(PK,Y2);
        PFC_LIN(i1,i2) = perfeval_pfc(PK,Y2);

        RUNTIMES(i1,i2) = toc;
    end
    runtime = ceil(sum(RUNTIMES(:))/60);
    minute_data = mod(runtime,60);
    hour_data = (runtime-minute_data)/60;
    fprintf(' %dh%dm\n',hour_data,minute_data);
end

clear SBANK0 SBANK1 SBANK2 Y1 Y2 U_MEAN U1 U2 R1 R2 W0 W1 PK

save(sprintf('%s_data.mat',mfilename));

figure('units','normalized','outerposition',[0 0 1 .6],'color','w')

subplot(1,2,1)
hold on
for i1 = 1:N1
    plot(thresholds,FC_LIN(i1,:),'color','b');
end
line([thresholds(1) thresholds(end)],[1/K 1/K],'color','k','linestyle',':');
axis([thresholds(1) thresholds(end) 0 1])
xlabel('threshold [mV]');
ylabel('performance');

subplot(1,2,2)
hold on
for i1 = 1:N1
    plot(thresholds,PFC_LIN(i1,:),'color','r');
end
line([thresholds(1) thresholds(end)],[1/K 1/K],'color','k','linestyle',':');
axis([thresholds(1) thresholds(end) 0 1])
xlabel('threshold [mV]');
ylabel('performance');

% show results:
saveas(gcf,sprintf('%s_done.png',mfilename));