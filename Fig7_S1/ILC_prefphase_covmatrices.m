clear all
close all

N = 500; M = 5000;

stim_DC = 0.5; V0 = -60; V1 = 12; noise_ref = 3;

lambda0 = 3; sigma0 = 2; c0 = 0.5;

cd(fileparts(which(mfilename)));
load('noise_levels.mat')

gx0 = {'const',0};
gy0 = {'const',0};
glambda = {'const',lambda0};
gtheta = {'grid',0,180,N};
gphi = {'const',0};
gsigma = {'const',sigma0};
GBANK = paramset(gx0,gy0,glambda,gtheta,gphi,gsigma);

s0 = [lambda0,90,0,stim_DC,stim_DC*c0];
SBANK0M = repmat(s0,M,1);
  
cov_matrix = NaN*zeros(N,N);

figure('units','normalized','outerposition',[0 0 1 1])

KMK = 100;
bin_width = 180/KMK;
bin_points = 50;
dx = bin_width/bin_points;

rows = 2;
P = length(thetasig_values);
for p = 1:P
    thetasig = thetasig_values(p);
    
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
    SBANK1 = SBANK0M;
    rands = rand(1,size(SBANK1,1));
    for sindex = 1:size(SBANK1,1)
        SBANK1(sindex,2) = SBANK1(sindex,2)+bin_indices(sum(ss < rands(sindex))+1)*bin_width;
    end
    
    noise = private_noise(p);
    U = V0+V1*cresponse(GBANK,SBANK1)+noise*randn(size(SBANK1,1),N);
    for i = 1:N
        for j = 1:N
            U1 = U(:,i);
            U2 = U(:,j);
            mean1 = mean(U1);
            mean2 = mean(U2);
            cov_matrix(i,j) = mean((U1-mean1).*(U2-mean2));
        end
    end
    subplot(rows,P/rows,p)
    imagesc(flipud(cov_matrix))
    daspect([1,1,1])
    box on
    set(gca,'xticklabel',[],'yticklabel',[])
    title({sprintf('thetasig = %g deg',thetasig),sprintf('private noise = %g mV',noise),''})
    colorbar('southoutside')
end

cd(fileparts(which(mfilename)));
saveas(gcf,sprintf('%s_N%d_refnoise%dmV_done.png',mfilename,N,noise_ref));