% script for generating data for Fig. 3
% generated data: performance curve of the linear decoder with phase nuisance
% extrenal functions used:
%   paramset
%   randisk
%   cresponse
%   mnrfitbb
%   mnrvals
%   nonlin
%   perfeval_fc
%   perfeval_pfc
%   maxfit

close all % close figures
clearvars % clear variables
cd(fileparts(which(mfilename))); % set current directory
addpath('../functions/'); % link external functions

N = 500; % number of cells (Gabor filters)
K = 10; % number of orientation categories
MK = 10; % number of stimulus orientation bins within a category (see M_theta in Table 1 in the paper)
J = 50; % number of bins for stimulus phase (see M_phi in Table 1 in the paper)
M1 = 20; % repetition of the stimulus bank for training (see M_rep in Table 1 in the paper)
M2 = 50; % repetition of the stimulus bank for independent test (see M_rep' in Table 1 in the paper)

% convergence settings for mnrfitbb (multinomial logistic regression)
maxiter = 5000; conv_crit = 1e-4; eta0 = 1e-8;

disk_radius = 90; % aperture radius in degrees (see R in Table 1 in the paper)
gcenter_seed = 1; % random seed for filter positions
gphi_seed = 1; % random seed for filter phases

stim_DC = 0.5;
V0 = -60; % mean membrane potential [mV] (see u_DC in Table 1 in the paper)
V1 = 12; % peak MP amplitude [mV] (see u_AC in Table 1 in the paper)
noise = 3; % std of Gaussian membrane potential noise [mV] (see sigma in Table 1 in the paper)

lambda0 = 3; % spatial period [deg]
sigma0 = 2; % size of the Gaussian envelope for circular Gabor filters [deg]
c0 = 0.5; % contrast

thresholds = -75:0.5:-40; %

N1 = length(thresholds);

% outputs:
FC_LIN = zeros(1,N1); % fraction correct
PFC_LIN = zeros(1,N1); % probablilistic fraction correct
RUNTIMES_LIN = zeros(1,N1); % run time of the simulation
ITERATIONS = zeros(1,N1); % number of iterations
% W_DATA = zeros(N1,N+1,K-1); % decoder weights if you want

% GABOR POPULATION:
  
gx0 = {'const',0}; % set later
gy0 = {'const',0}; % set later
glambda = {'const',lambda0}; % all filters has the same spatial period
gtheta = {'grid',0,180,N}; % evenly distributed orientations
gphi = {'perm',0,360,gphi_seed}; % evenly distributed phases but order is shuffled
gsigma = {'const',sigma0}; % all filters has the same Gaussian envelope
GBANK = paramset(gx0,gy0,glambda,gtheta,gphi,gsigma);
GBANK(:,1:2) = randisk(disk_radius,N,gcenter_seed); % filter positions: uniform distribution on a disk

% STIMULUS BANK:

s0lambda = {'const',lambda0}; % fix spatial period and the same as glambda
s0theta = {'grid',0,180,K*MK}; % stimulus orientiation grid
s0phi = {'grid',0,360,J}; % stimulus phase grid
s0DC = {'const',stim_DC}; % fix DC component
s0AC = {'const',stim_DC*c0}; % fix contrast

SBANK0 = paramset(s0lambda,s0theta,s0phi,s0DC,s0AC);
SBANK1 = repmat(SBANK0,M1,1); % for training
SBANK2 = repmat(SBANK0,M2,1); % for testing

% ORIENTATION CATEGORIES:

Y1 = ceil(SBANK1(:,2)/(180/K)); % orienation lables: 1,2,...,K
Y2 = ceil(SBANK2(:,2)/(180/K)); % orienation lables: 1,2,...,K

% MP RESPONSE:

U_MEAN = V0+V1*cresponse(GBANK,SBANK0); % filter response

% adding noise
U1 = repmat(U_MEAN,M1,1)+noise*randn(K*MK*J*M1,N); % for training
U2 = repmat(U_MEAN,M2,1)+noise*randn(K*MK*J*M2,N); % for testing

W0 = 0; % initial weights
for i1 = 1:N1 % go through FRNL thresholds

    threshold = thresholds(i1); % FRNL threshold
    
    tic % starting timer
    
    R1 = nonlin(U1,threshold,1); % firing rates for training
    R2 = nonlin(U2,threshold,1); % firing rates for testing
    
    [W1,etas] = mnrfitbb(R1,Y1,W0,maxiter,conv_crit,eta0); % find the best linear weights with gradient descent
    ITERATIONS(i1) = length(etas); % number of iterations used
    % W_DATA(i1,:,:) = W1; % save decoder weights if you want
    W0 = W1; % setting initial condition for the next threshold
    
    if i1 == 1 fprintf('%d',ITERATIONS(i1)); else fprintf('-%d',ITERATIONS(i1)); end % log print
    
    PK_LIN = mnrvals(W1,R2); % posterior probabilities
    
    FC_LIN(i1) = perfeval_fc(PK_LIN,Y2); % calculating fraction correct 
    PFC_LIN(i1) = perfeval_pfc(PK_LIN,Y2); % calculating probabilistic fraction correct
    
    RUNTIMES_LIN(i1) = toc; % stop timer
end

runtime = ceil(sum(RUNTIMES_LIN(:))/60); % sum all times
minute_data = mod(runtime,60);
hour_data = (runtime-minute_data)/60;
fprintf(' %dh%dm\n',hour_data,minute_data);

clear SBANK0 SBANK1 SBANK2 Y1 Y2 U_MEAN U1 U2 R1 R2 W0 W1 PK_LIN % clear not important large data

save(sprintf('%s_data.mat',mfilename)); % save all remaining important data

% make a figure just for check
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
plot([])
hold on
axis([thresholds(1) thresholds(end) 0 1])
xlabel('threshold [mV]');
ylabel('performance');
title({sprintf('N = %d, K = %d, MK = %d, J = %d, M1 = %d, M2 = %d',N,K,MK,J,M1,M2),sprintf('convcrit = %g, eta0 = %g, resolution = %g [mV] (runtime: %d)',conv_crit,eta0,thresholds(2)-thresholds(1),runtime) })
p1 = plot(thresholds,FC_LIN,'color','b','linewidth',2);
p2 = plot(thresholds,PFC_LIN,'color','b','linewidth',2,'linestyle',':');
l0 = line([thresholds(1) thresholds(end)],[1/K 1/K],'color','r','linestyle',':');
[maxx_fc,maxy_fc] = maxfit(thresholds,FC_LIN,3);
[maxx_pfc,maxy_pfc] = maxfit(thresholds,PFC_LIN,3);
line([maxx_fc maxx_fc],[0 0.033],'color','b','linewidth',2)
line([maxx_pfc maxx_pfc],[0 0.033],'color','b','linewidth',2,'linestyle',':')
scatter(maxx_fc,maxy_fc,100,'r','d','filled')
scatter(maxx_pfc,maxy_pfc,100,'r','d','filled')
legend([p1,p2,l0],{'FC','PFC','chance level'})

saveas(gcf,sprintf('%s_done.png',mfilename)); % save figure