% script for generating data for Fig. 3
% generated data: performance curve for pixel and DoG populations (with phase nuisance)
% extrenal functions used:
%   paramset
%   randisk
%   cresponse
%   mnrfitbb
%   mnrvals
%   nonlin
%   perfeval_fc
%   perfeval_pfc

% see Table 1 in the paper (MK: M_theta, disk_radius: R, V0: u_DC, V1: u_AC)
N = 500; K = 10; MK = 10; J = 50; M1 = 10; M2 = 20;

maxiter = 5000; conv_crit = 1e-4; eta0 = 1e-8;

lambda0 = 3; stim_DC = 0.5; c0 = 0.5;

disk_radius = 90; center_seed = 1;

Gabor_lambda = lambda0; Gabor_sigma = 2; gphi_seed = 1;

DoG_sigma1 = 1.5; DoG_sigma2 = 4; DoG_lambda = 1000;

V0 = -60; V1 = 12; noise = 3; 

thresholds = -75:0.5:-40;

N1 = length(thresholds);

FC_Gabor = zeros(1,N1);
PFC_Gabor = zeros(1,N1);
RUNTIMES_Gabor = zeros(1,N1);
ITERATIONS_Gabor = zeros(1,N1);
WDATA_Gabor = zeros(N1,N+1,K-1);

FC_pixel = zeros(1,N1);
PFC_pixel = zeros(1,N1);
RUNTIMES_pixel = zeros(1,N1);
ITERATIONS_pixel = zeros(1,N1);
WDATA_pixel = zeros(N1,N+1,K-1);

FC_DoG = zeros(1,N1);
PFC_DoG = zeros(1,N1);
RUNTIMES_DoG = zeros(1,N1);
ITERATIONS_DoG = zeros(1,N1);
WDATA_DoG = zeros(N1,N+1,K-1);

% GABOR POPULATION:
  
gx0 = {'const',0};
gy0 = {'const',0};
glambda = {'const',Gabor_lambda};
gtheta = {'grid',0,180,N};
gphi = {'perm',0,360,gphi_seed};
gsigma = {'const',Gabor_sigma};
GBANK = paramset(gx0,gy0,glambda,gtheta,gphi,gsigma);
GBANK(:,1:2) = randisk(disk_radius,N,center_seed);

X0 = GBANK(:,1)';
Y0 = GBANK(:,2)';

% DOG POPULATION:

GBANK1 = GBANK;
GBANK2 = GBANK;
GBANK1(:,3) = DoG_lambda;
GBANK2(:,3) = DoG_lambda;
GBANK1(:,4:5) = 0;
GBANK2(:,4:5) = 0;
GBANK1(:,6) = DoG_sigma1;
GBANK2(:,6) = DoG_sigma2;

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

UMEAN_Gabor = V0+V1*cresponse(GBANK,SBANK0);
Gabor_std = std(UMEAN_Gabor(:));

UMEAN_DoG_unscaled = cresponse(GBANK1,SBANK0)-cresponse(GBANK2,SBANK0); % DoG response is calculated as difference of to degenerated Gabor responses
UMEAN_pixel_unscaled = (1+cos((2*pi*(cos(pi*SBANK0(:,2)/180)*Y0-sin(pi*SBANK0(:,2)/180)*X0)./repmat(SBANK0(:,1),1,N))-pi*SBANK0(:,3)/180))/2;

DoG_meanvec2d = mean(UMEAN_DoG_unscaled);
DoG_stdvec2d =  std(UMEAN_DoG_unscaled);
pixel_meanvec2d = mean(UMEAN_pixel_unscaled);
pixel_stdvec2d =  std(UMEAN_pixel_unscaled);

M0 = size(SBANK0,1);

UMEAN_DoG_scaled = V0+Gabor_std*(UMEAN_DoG_unscaled-repmat(DoG_meanvec2d,M0,1))./repmat(DoG_stdvec2d,M0,1);
UMEAN_pixel_scaled = V0+Gabor_std*(UMEAN_pixel_unscaled-repmat(pixel_meanvec2d,M0,1))./repmat(pixel_stdvec2d,M0,1);

U1_Gabor = repmat(UMEAN_Gabor,M1,1)+noise*randn(M0*M1,N);
U2_Gabor = repmat(UMEAN_Gabor,M2,1)+noise*randn(M0*M2,N);

U1_DoG = repmat(UMEAN_DoG_scaled,M1,1)+noise*randn(M0*M1,N);
U2_DoG = repmat(UMEAN_DoG_scaled,M2,1)+noise*randn(M0*M2,N);

U1_pixel = repmat(UMEAN_pixel_scaled,M1,1)+noise*randn(M0*M1,N);
U2_pixel = repmat(UMEAN_pixel_scaled,M2,1)+noise*randn(M0*M2,N);

W0_Gabor = 0;
W0_DoG = 0;
W0_pixel = 0;
for i1 = 1:N1

    threshold = thresholds(i1);
    
    R1_Gabor = nonlin(U1_Gabor,threshold,1);
    R2_Gabor = nonlin(U2_Gabor,threshold,1);
    
    R1_DoG = nonlin(U1_DoG,threshold,1);
    R2_DoG = nonlin(U2_DoG,threshold,1);
    
    R1_pixel = nonlin(U1_pixel,threshold,1);
    R2_pixel = nonlin(U2_pixel,threshold,1);
    
    tic
    
    [W1_Gabor,etas] = mnrfitbb(R1_Gabor,Y1,W0_Gabor,maxiter,conv_crit,eta0);
    ITERATIONS_Gabor(i1) = length(etas);
    WDATA_Gabor(i1,:,:) = W1_Gabor;
    W0_Gabor = W1_Gabor;
    RUNTIMES_Gabor(i1) = toc;
    if i1 == 1 fprintf('%d',ITERATIONS_Gabor(i1)); else fprintf('-%d',ITERATIONS_Gabor(i1)); end
    PK = mnrvals(W1_Gabor,R2_Gabor);
    FC_Gabor(i1) = perfeval_fc(PK,Y2);
    PFC_Gabor(i1) = perfeval_pfc(PK,Y2);
    
    tic
    
    [W1_DoG,etas] = mnrfitbb(R1_DoG,Y1,W0_DoG,maxiter,conv_crit,eta0);
    ITERATIONS_DoG(i1) = length(etas);
    WDATA_DoG(i1,:,:) = W1_DoG;
    W0_DoG = W1_DoG;
    RUNTIMES_DoG(i1) = toc;
    fprintf(',%d',ITERATIONS_DoG(i1))
    PK = mnrvals(W1_DoG,R2_DoG);
    FC_DoG(i1) = perfeval_fc(PK,Y2);
    PFC_DoG(i1) = perfeval_pfc(PK,Y2);
    
    tic
    
    [W1_pixel,etas] = mnrfitbb(R1_pixel,Y1,W0_pixel,maxiter,conv_crit,eta0);
    ITERATIONS_pixel(i1) = length(etas);
    WDATA_pixel(i1,:,:) = W1_pixel;
    W0_pixel = W1_pixel;
    RUNTIMES_pixel(i1) = toc;
    fprintf(',%d',ITERATIONS_pixel(i1))
    PK = mnrvals(W1_pixel,R2_pixel);
    FC_pixel(i1) = perfeval_fc(PK,Y2);
    PFC_pixel(i1) = perfeval_pfc(PK,Y2);
end
fprintf('\n');

clear SBANK0 SBANK1 SBANK2 Y1 Y2 UMEAN_Gabor UMEAN_DoG_unscaled UMEAN_pixel_unscaled UMEAN_DoG_scaled UMEAN_pixel_scaled U1_Gabor U2_Gabor U1_DoG U2_DoG U1_pixel U2_pixel R1_Gabor R2_Gabor R1_DoG R2_DoG R1_pixel R2_pixel W0_Gabor W0_DoG W0_pixel W1_Gabor W1_DoG W1_pixel PK

save(sprintf('%s_data.mat',mfilename));

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
plot([])
hold on
axis([thresholds(1) thresholds(end) 0 1])
xlabel('threshold [mV]');
ylabel('performance');
title({sprintf('N = %d, K = %d, MK = %d, J = %d, M1 = %d, M2 = %d',N,K,MK,J,M1,M2),sprintf('convcrit = %g, eta0 = %g, resolution = %g mV',conv_crit,eta0,thresholds(2)-thresholds(1)) })
p1 = plot(thresholds,FC_Gabor,'color','b','linewidth',2);
p2 = plot(thresholds,FC_DoG,'color',[0 .5 0],'linewidth',2);
p3 = plot(thresholds,FC_pixel,'color','r','linewidth',2);
p4 = plot(thresholds,PFC_Gabor,'color','b','linewidth',2,'linestyle',':');
p5 = plot(thresholds,PFC_DoG,'color',[0 .5 0],'linewidth',2,'linestyle',':');
p6 = plot(thresholds,PFC_pixel,'color','r','linewidth',2,'linestyle',':');
l0 = line([thresholds(1) thresholds(end)],[1/K 1/K],'color','r','linestyle',':');
legend([p1,p4,p2,p5,p3,p6,l0],{'Gabor FC','Gabor PFC','DoG FC','DoG PFC','pixel FC','pixel PFC','chance level'})

saveas(gcf,sprintf('%s_done.png',mfilename));