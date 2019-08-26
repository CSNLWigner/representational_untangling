% script for generating data for Fig. 4 figure supplement 2

close all
clearvars

N = 500; K = 10; MK = 10; J = 50;

M1 = 100000; M2 = 200000;

maxiter = 10000; conv_crit = 1e-3; eta0 = 1e-8;

nuisance(1).name = 'phase';
nuisance(1).label = '$$\varphi$$';
nuisance(1).range = [0 360];
nuisance(1).mean = 180;

nuisance(2).name = 'wavelength';
nuisance(2).label = '$$\lambda$$';
nuisance(2).range = [0.5 7.5];
nuisance(2).mean = 3;

nuisance(3).name = 'contrast';
nuisance(3).label = '$$c$$';
nuisance(3).range = [0 1];
nuisance(3).mean = 0.4;

phi0 = 180;
lambda0 = 3;
c0 = 0.5; 

pmu = 0.95;
psigma = 0.55;
lognorm = @(x) exp(-(log(x)-pmu).^2/psigma^2/2)./x/psigma/sqrt(2*pi);

palpha = 2.4;
pbeta = 3.6;
beta = @(x) (x.^(palpha-1)).*(1-x).^(pbeta-1);

boundaries_wavelength = bdd_splitter(nuisance(2).range(1),nuisance(2).range(2),lognorm,K);
boundaries_contrast = bdd_splitter(nuisance(3).range(1),nuisance(3).range(2),beta,K);

disk_radius = 90; sigma0 = 2; stim_DC = 0.5; V0 = -60; V1 = 12; noise = 3; 

gseeds = 1:25;

thresholds = -60:1:-50;

N1 = length(gseeds);
N2 = length(thresholds);

FC = zeros(N1,N2);
ITERATIONS = zeros(N1,N2);
RUNTIMES = zeros(N1,N2);

load('SBANK.mat','SBANK1','SBANK2')

S1 = repmat(SBANK1(:,5),1,K)/stim_DC;
S2 = repmat(SBANK2(:,5),1,K)/stim_DC;
B1 = repmat(boundaries_contrast(1:end-1),M1,1);
B2 = repmat(boundaries_contrast(1:end-1),M2,1);
Y1 = sum((S1 > B1)')';
Y2 = sum((S2 > B2)')';

warning off

for i1 = 1:N1

gx0 = {'const',0};
gy0 = {'const',0};
glambda = {'const',lambda0};
gtheta = {'grid',0,180,N};
gphi = {'perm',0,360,gseeds(i1)};
gsigma = {'const',sigma0};
GBANK = paramset(gx0,gy0,glambda,gtheta,gphi,gsigma);
GBANK(:,1:2) = randisk(disk_radius,N,gseeds(i1));
GBANK(:,3) = samples(nuisance(2).range(1),nuisance(2).range(2),lognorm,N,'permutated_representative',gseeds(i1));

U1 = V0+V1*cresponse(GBANK,SBANK1)+noise*randn(M1,N);
U2 = V0+V1*cresponse(GBANK,SBANK2)+noise*randn(M2,N);

W0 = 0;
for i2 = 1:N2
    
    threshold = thresholds(i2);
    
    tic;
    
    R1 = nonlin(U1,threshold,1);
    R2 = nonlin(U2,threshold,1);
    
    [W1,etas] = mnrfitbb(R1,Y1,W0,maxiter,conv_crit,eta0);
    ITERATIONS(i1,i2) = length(etas);
    W0 = W1;
    
    if i2 == 1 fprintf('%d',ITERATIONS(i1,i2)); else fprintf('-%d',ITERATIONS(i1,i2)); end
    
    PK = mnrvals(W1,R2);
    
    FC(i1,i2) = perfeval_fc(PK,Y2);
    
    RUNTIMES(i1,i2) = toc/60;
end

runtime = ceil(sum(RUNTIMES(i1,:))/60);
minute_data = mod(runtime,60);
hour_data = (runtime-minute_data)/60;
fprintf(' /%dh%dm/\n',hour_data,minute_data);

end

cd(fileparts(which(mfilename)));
clear SBANK1 SBANK2 S1 S2 B1 B2 Y1 Y2 U1 U2 R1 R2 PK

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
    curve_plots(i1) = plot(thresholds,FC(i1,:),'color',color);
    legend_strings{end+1} = sprintf('gseed = %d (%dmin)',gseeds(i1),ceil(sum(RUNTIMES(i1,:))/60));
    [~,maxi] = max(FC(i1,:)); maxx0 = thresholds(maxi);
    di = 3; di1 = di; di2 = di;
    params = polyfit(thresholds(maxi-di1:maxi+di2),FC(i1,maxi-di1:maxi+di2),2);
    maxx = -params(2)/params(1)/2;
    if abs(abs(maxx-maxx0)-dth/2) < dth/6
        di1 = di-(sign(maxx-maxx0)+1)/2;
        di2 = di-1+(sign(maxx-maxx0)+1)/2;
        params = polyfit(thresholds(maxi-di1:maxi+di2),FC(i1,maxi-di1:maxi+di2),2);
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

save([mfilename,'_data.mat']);
saveas(gcf,sprintf('%s_done.png',mfilename));