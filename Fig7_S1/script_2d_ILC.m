clear all
close all

N = 500; K = 10; MK = 10; J = 50; M1 = 10; M2 = 20;

maxiter = 5000; conv_crit = 1e-4; eta0 = 1e-8;

disk_radius = 90; gcenter_seed = 1; gphi_seed = 1;

stim_DC = 0.5; V0 = -60; V1 = 12; 

cd(fileparts(which(mfilename)));
load('noise_levels.mat')

thetasig_index = 5;
thetasig = thetasig_values(thetasig_index);
noise = private_noise(thetasig_index); 

lambda0 = 3; sigma0 = 2; c0 = 0.5;

thresholds = -75:0.5:-40;

N1 = length(thresholds);

FC_LIN = zeros(1,N1);
PFC_LIN = zeros(1,N1);
RUNTIMES_LIN = zeros(1,N1);
ITERATIONS = zeros(1,N1);

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

bin_width = 180/K/MK;
bin_points = 50;
dx = bin_width/bin_points;
bins = round(10*thetasig/bin_width);
if mod(bins,2) == 0, bins = bins+1; end
xmin = -bins*bin_width/2+dx/2;
xmax = bins*bin_width/2-dx/2;
xx = xmin:dx:xmax+dx/2;
yy = exp(-xx.^2/thetasig^2/2);
yy = yy/sum(yy);
bin_indices = -(bins-1)/2:(bins-1)/2;
probabilities = zeros(1,bins);
for k = 1:bins
    probabilities(k) = sum(yy((k-1)*bin_points+1:(k-1)*bin_points+bin_points));
end
ss = cumsum(probabilities);

rands = rand(1,size(SBANK1,1));
for sindex = 1:size(SBANK1,1)
    SBANK1(sindex,2) = SBANK1(sindex,2)+bin_indices(sum(ss < rands(sindex))+1)*bin_width;
end
rands = rand(1,size(SBANK2,1));
for sindex = 1:size(SBANK2,1)
    SBANK2(sindex,2) = SBANK2(sindex,2)+bin_indices(sum(ss < rands(sindex))+1)*bin_width;
end

U1 = V0+V1*cresponse(GBANK,SBANK1)+noise*randn(size(SBANK1,1),N);
U2 = V0+V1*cresponse(GBANK,SBANK2)+noise*randn(size(SBANK2,1),N);

[counts,centers] = hist(U1(:),100);
save(sprintf('Uhist_data.mat',mfilename),'counts','centers');

for i1 = 1:N1

    threshold = thresholds(i1);
    
    tic
    
    R1 = nonlin(U1,threshold,1);
    R2 = nonlin(U2,threshold,1);
    
    [W1,etas] = mnrfitbb(R1,Y1,0,maxiter,conv_crit,eta0);
    ITERATIONS(i1) = length(etas);
    
    if i1 == 1 fprintf('%d',ITERATIONS(i1)); else fprintf('-%d',ITERATIONS(i1)); end
    
    PK_LIN = mnrvals(W1,R2);
    
    FC_LIN(i1) = perfeval_fc(PK_LIN,Y2);
    PFC_LIN(i1) = perfeval_pfc(PK_LIN,Y2);
    
    RUNTIMES_LIN(i1) = toc;
end

runtime = ceil(sum(RUNTIMES_LIN(:))/60);
minute_data = mod(runtime,60);
hour_data = (runtime-minute_data)/60;
fprintf(' %dh%dm\n',hour_data,minute_data);

clear SBANK0 SBANK1 SBANK2 Y1 Y2 U1 U2 R1 R2 W1 PK_LIN

save(sprintf('%s_data.mat',mfilename));

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

saveas(gcf,sprintf('%s_done.png',mfilename));