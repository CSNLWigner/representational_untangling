% make 'noiselevels.mat' and other data in folder 'script_noise_scans'
close all
clearvars

cd(fileparts(which(mfilename)));

N = 500; K = 10; MK = 10; J = 50; MS = 50;

disk_radius = 90; gcenter_seed = 1; gphi_seed = 1;

stim_DC = 0.5; V0 = -60; V1 = 12;

lambda0 = 3; sigma0 = 2; c0 = 0.5;

thetasig_values = 1:12;
N1 = length(thetasig_values);

noise_min = 2.55;
noise_max = 2.80;
noise_step = 0.01;
noise_values = noise_min:noise_step:noise_max;
N2 = length(noise_values);

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

bin_width = 180/K/MK;
bin_points = 50;
dx = bin_width/bin_points;
private_noise = zeros(1,N1);
for n1 = 1:N1
    thetasig = thetasig_values(n1);
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
    SBANK1 = repmat(SBANK0,MS,1);
    rands = rand(1,size(SBANK1,1));
    for sindex = 1:size(SBANK1,1)
        SBANK1(sindex,2) = SBANK1(sindex,2)+bin_indices(sum(ss < rands(sindex))+1)*bin_width;
    end
    
    VARS = zeros(1,N2);
    for n2 = 1:N2
        noise = noise_values(n2);
        U = V0+V1*cresponse(GBANK,SBANK1)+noise*randn(size(SBANK1,1),N);
        U3D = reshape(U',N,size(SBANK0,1),MS);
        VARS(n2) = mean(mean(var(U3D,0,3)));
    end

    y0 = 9;

    i1 = sum(VARS < y0);
    i2 = i1+1;

    x1 = noise_values(i1);
    x2 = noise_values(i2);
    y1 = VARS(i1);
    y2 = VARS(i2);

    a = (y2-y0)/(y2-y1);
    x = a*x1+(1-a)*x2;

    private_noise(n1) = x;
    
    hold on
    plot(noise_values,VARS,'color','b')
    plot(noise_values(i1:i2),VARS(i1:i2),'color','b','linewidth',2)
    scatter(noise_values,VARS,30,'b','filled')
    line([noise_values(1) noise_values(end)],[y0 y0],'color','k')

    line([x x],[VARS(1) VARS(end)],'color','r','linewidth',2)
    xlim([noise_values(1) noise_values(end)])
    ylim([VARS(1) VARS(end)])
    xlabel('private noise [mV]')
    ylabel('mean combined variance [mV^2]')
    title({sprintf('thetasig = %d deg, red solution at %g mV',thetasig,x),''})

    saveas(gcf,sprintf('script_noise_scans/%s_thetasig%d.png',mfilename,thetasig));
    close all
end
save('noise_levels.mat','thetasig_values','private_noise')