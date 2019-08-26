% script for generating data for Fig. 3
% generated data: sparsness data (with phase nuisance)
% extrenal functions used:
%   paramset
%   randisk
%   cresponse
%   nonlin

% see Table 1 in the paper (MK: M_theta, disk_radius: R, V0: u_DC, V1: u_AC)
N = 500; K = 10; MK = 10; J = 50; M = 100;

disk_radius = 90; gcenter_seed = 1; gphi_seed = 1;

stim_DC = 0.5; V0 = -60; V1 = 12; noise = 3; 

lambda0 = 3; sigma0 = 2; c0 = 0.5;

thresholds = -75:0.5:-40;

N1 = length(thresholds);

SPARSENESS = zeros(1,N1);

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
SBANK = repmat(SBANK0,M,1);

KMKJM = K*MK*J*M;

U = V0+V1*cresponse(GBANK,SBANK)+noise*randn(KMKJM,N);

RHAT_NAN = [];
SPPOP_NAN = [];
for i1 = 1:N1

    threshold = thresholds(i1);
    
    R = nonlin(U,threshold,1);
    
    RHAT = R./repmat(sqrt(sum(R.^2)/KMKJM),KMKJM,1); % Eq. S35
    
    SP_POP = (1-(sum(RHAT')/N).^2./(sum(RHAT'.^2)/N))/(1-1/N); % Eq. S36
    
    RHAT_NAN(end+1) = sum(sum(isnan(RHAT)));
    SPPOP_NAN(end+1) = sum(isnan(SP_POP));
    
    SPARSENESS(i1) = sum(SP_POP(~isnan(SP_POP)))/length(SP_POP(~isnan(SP_POP)));
    
    fprintf('.');
end

fprintf('\n');

clear SBANK0 SBANK U R

save(sprintf('%s_data.mat',mfilename));

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
plot([])
hold on
axis([thresholds(1) thresholds(end) 0 1])
xlabel('threshold [mV]');
ylabel('sparseness');
title(sprintf('N = %d, K = %d, MK = %d, J = %d, M = %d',N,K,MK,J,M))
sparseness_color = [255 135 137]/255;
plot(thresholds,SPARSENESS,'color',sparseness_color,'linewidth',3)

saveas(gcf,sprintf('%s_done.png',mfilename));