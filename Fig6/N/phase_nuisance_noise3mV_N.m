clear all
close all

N_values = [125,250,500,1000];

N_reference = 500; reference_index = find(N_values == N_reference);

K = 10; MK = 10; J = 50; M1 = 20; M2 = 50;

maxiter = 10000; conv_crit = 1e-4; eta0 = 1e-8;

disk_radius = 90; gcenter_seed = 1; gphi_seed = 1;

stim_DC = 0.5; V0 = -60; V1 = 12; noise = 3; 

lambda0 = 3; sigma0 = 2; c0 = 0.5;

thresholds = -75:0.5:-40;

N1 = length(N_values);
N2 = length(thresholds);

FC_LIN = zeros(N1,N2);
PFC_LIN = zeros(N1,N2);
RUNTIMES_LIN = zeros(N1,N2);
ITERATIONS = zeros(N1,N2);

for i1 = 1:N1
    
    fprintf('curve %d/%d:\n',i1,N1);
    
    N = N_values(i1);
    
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
    
    PK_LIN = mnrvals(W1,R2);
    
    FC_LIN(i1,i2) = perfeval_fc(PK_LIN,Y2);
    PFC_LIN(i1,i2) = perfeval_pfc(PK_LIN,Y2);
    
    RUNTIMES_LIN(i1,i2) = toc;
end

runtime = ceil(sum(RUNTIMES_LIN(i1,:))/60);
minute_data = mod(runtime,60);
hour_data = (runtime-minute_data)/60;
fprintf(' /%dh%dm/\n',hour_data,minute_data);

end

cd(fileparts(which(mfilename)));
save(sprintf('%s_data.mat',mfilename),'N_values','N_reference','reference_index','K','MK','J','M1','M2','maxiter','conv_crit','eta0','disk_radius','gcenter_seed','gphi_seed','stim_DC','V0','V1','noise','lambda0','sigma0','c0','thresholds','FC_LIN','PFC_LIN','RUNTIMES_LIN','ITERATIONS','GBANK');

MAXX = zeros(1,N1);
MAXY = zeros(1,N1);

figure('units','normalized','outerposition',[0 0 1 1])
plot([])
hold on
axis([thresholds(1) thresholds(end) 0 1])
line([thresholds(1) thresholds(end)],[1/K 1/K],'linestyle',':','color','r')
xlabel('threshold [mV]');
ylabel('probabilistic fraction correct');
dth = thresholds(2)-thresholds(1);
title({sprintf('K = %d, MK = %d, J = %d, M1 = %d, M2 = %d',K,MK,J,M1,M2),sprintf('convcrit = %g, eta0 = %g, resolution = %g [mV]',conv_crit,eta0,dth)})
curve_plots = 1:N1;
legend_strings = {};
for i1 = 1:N1
    if i1 == reference_index, ref = ' *'; else, ref = ''; end
    color = [0 1-i1/N1 i1/N1];
    curve_plots(i1) = plot(thresholds,PFC_LIN(i1,:),'color',color);
    legend_strings{end+1} = sprintf('N = %d (%dmin)%s',N_values(i1),ceil(sum(RUNTIMES_LIN(i1,:))/60),ref);
    [~,maxi] = max(PFC_LIN(i1,:)); maxx0 = thresholds(maxi);
    di = 3; di1 = di; di2 = di;
    params = polyfit(thresholds(maxi-di1:maxi+di2),PFC_LIN(i1,maxi-di1:maxi+di2),2);
    maxx = -params(2)/params(1)/2;
    if abs(abs(maxx-maxx0)-dth/2) < dth/6
        di1 = di-(sign(maxx-maxx0)+1)/2;
        di2 = di-1+(sign(maxx-maxx0)+1)/2;
        params = polyfit(thresholds(maxi-di1:maxi+di2),PFC_LIN(i1,maxi-di1:maxi+di2),2);
        maxx = -params(2)/params(1)/2;
    end
    maxy = params(3)-params(2)^2/params(1)/4;
    scatter(maxx,maxy,10,color,'filled')
    line([maxx maxx],[0 0.033],'color',color)
    xfit = linspace(thresholds(maxi-di1),thresholds(maxi+di2),100);
    yfit = params(1)*xfit.^2+params(2)*xfit+params(3);
    plot(xfit,yfit,'color','r');
    MAXX(i1) = maxx;
    MAXY(i1) = maxy;
end
legend(curve_plots,legend_strings)
legend boxoff

save(sprintf('%s_data.mat',mfilename),'MAXX','MAXY','-append');

saveas(gcf,sprintf('%s_done.png',mfilename));