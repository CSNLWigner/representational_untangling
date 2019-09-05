close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../../functions/');

kappa = 1;

noise_values = 1.5:0.5:7;

N = 500; K = 10; MK = 10; J = 50; M1 = 25; M2 = 50;

maxiter = 10000; conv_crit = 1e-4; eta0 = 1e-8;

disk_radius = 90; gcenter_seed = 1; gphi_seed = 1;

stim_DC = 0.5; V0 = -60; V1 = 12;

lambda0 = 3; sigma0 = 2; c0 = 0.5;

thresholds = -65:0.5:-45;

N1 = length(noise_values);
N2 = length(thresholds);

FC = zeros(N1,N2);
PFC = zeros(N1,N2);
RUNTIMES = zeros(N1,N2);
ITERATIONS = zeros(N1,N2);
    
% GABOR POPULATION:

gx0 = {'const',0};
gy0 = {'const',0};
glambda = {'const',lambda0};
gtheta = {'grid',0,180,N};
gphi = {'perm',0,360,gphi_seed};
gsigma = {'const',sigma0};
GBANK = paramset(gx0,gy0,glambda,gtheta,gphi,gsigma);
GBANK(:,1:2) = randisk(disk_radius,N,gcenter_seed);

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

for i1 = 1:N1
    
    fprintf('curve %d/%d:\n',i1,N1);
    
    noise = noise_values(i1);
    
    U1 = repmat(U_MEAN,M1,1)+noise*randn(K*MK*J*M1,N);
    U2 = repmat(U_MEAN,M2,1)+noise*randn(K*MK*J*M2,N);
    
for i2 = 1:N2
        
    threshold = thresholds(i2);
    
    tic
    
    R1 = nonlin(U1,threshold,1,kappa);
    R2 = nonlin(U2,threshold,1,kappa);
    
    [W1,etas] = mnrfitbb(R1,Y1,0,maxiter,conv_crit,eta0);
    ITERATIONS(i1,i2) = length(etas);
    
    if i2 == 1 fprintf('%d',ITERATIONS(i1,i2)); else fprintf('-%d',ITERATIONS(i1,i2)); end
    
    PK = mnrvals(W1,R2);
    
    FC(i1,i2) = perfeval_fc(PK,Y2);
    PFC(i1,i2) = perfeval_pfc(PK,Y2);
    
    RUNTIMES(i1,i2) = toc;
end

curve_runtime = ceil(sum(RUNTIMES(i1,:))/60);
minute_data = mod(curve_runtime,60);
hour_data = (curve_runtime-minute_data)/60;
fprintf(' /%dh%dm/\n',hour_data,minute_data);

end

save(sprintf('%s_data.mat',mfilename),'kappa','noise_values','N','K','MK','J','M1','M2','maxiter','conv_crit','eta0','disk_radius','gcenter_seed','gphi_seed','stim_DC','V0','V1','lambda0','sigma0','c0','thresholds','FC','PFC','RUNTIMES','ITERATIONS','GBANK');

MAXX = zeros(1,N1);
MAXY = zeros(1,N1);

figure('units','normalized','outerposition',[0 0 1 1])
plot([])
hold on
axis([thresholds(1) thresholds(end) 0 1])
line([thresholds(1) thresholds(end)],[1/K 1/K],'linestyle',':','color','r')
xlabel('threshold [mV]');
ylabel('fraction correct');
dth = thresholds(2)-thresholds(1);
title({sprintf('N = %d, K = %d, MK = %d, J = %d, M1 = %d, M2 = %d',N,K,MK,J,M1,M2),sprintf('convcrit = %g, eta0 = %g, resolution = %g [mV]',conv_crit,eta0,dth)})
curve_plots = 1:N1;
legend_strings = {};
for i1 = 1:N1
    color = [0 1-i1/N1 i1/N1];
    [maxx,maxy,xfit,yfit] = maxfit(thresholds,FC(i1,:),3,50);
    line([maxx maxx],[0 0.033],'color',color)
    scatter(maxx,maxy,10,color,'filled')
    plot(xfit,yfit,'color','r','linewidth',2);
    MAXX(i1) = maxx;
    MAXY(i1) = maxy; 
    curve_plots(i1) = plot(thresholds,FC(i1,:),'color',color,'linewidth',1);
    legend_strings{end+1} = sprintf('noise = %.2g (%dmin)',noise_values(i1),ceil(sum(RUNTIMES(i1,:))/60));
end
legend(curve_plots,legend_strings)
legend boxoff

save(sprintf('%s_data.mat',mfilename),'MAXX','MAXY','-append');

saveas(gcf,sprintf('%s_done.png',mfilename));