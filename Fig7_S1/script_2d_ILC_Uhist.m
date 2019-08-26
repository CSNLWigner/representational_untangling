clear all
close all

N = 500; K = 10; MK = 10; J = 50; M1 = 25;

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

U1 = V0+V1*cresponse(GBANK,SBANK1)+noise*randn(size(SBANK1,1),N);

[counts,centers] = hist(U1(:),100);

save(sprintf('Uhist_data.mat',mfilename),'counts','centers');