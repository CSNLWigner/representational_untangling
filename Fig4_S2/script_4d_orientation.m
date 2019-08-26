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

% boundaries_wavelength = bdd_splitter(nuisance(2).range(1),nuisance(2).range(2),lognorm,K);
% boundaries_contrast = bdd_splitter(nuisance(3).range(1),nuisance(3).range(2),beta,K);

disk_radius = 90; sigma0 = 2; stim_DC = 0.5; V0 = -60; V1 = 12; noise = 3; 

gseed = 25;

thresholds = -75:.5:-40;

N2 = length(thresholds);

FC = zeros(1,N2);
ITERATIONS = zeros(1,N2);
RUNTIMES = zeros(1,N2);

load('SBANK.mat','SBANK1','SBANK2')

Y1 = ceil(SBANK1(:,2)/(180/K));
Y2 = ceil(SBANK2(:,2)/(180/K));

gx0 = {'const',0};
gy0 = {'const',0};
glambda = {'const',lambda0};
gtheta = {'grid',0,180,N};
gphi = {'perm',0,360,gseed};
gsigma = {'const',sigma0};
GBANK = paramset(gx0,gy0,glambda,gtheta,gphi,gsigma);
GBANK(:,1:2) = randisk(disk_radius,N,gseed);
GBANK(:,3) = samples(nuisance(2).range(1),nuisance(2).range(2),lognorm,N,'permutated_representative',gseed);

U1 = V0+V1*cresponse(GBANK,SBANK1)+noise*randn(M1,N);
U2 = V0+V1*cresponse(GBANK,SBANK2)+noise*randn(M2,N);

warning off

W0 = 0;
for i2 = 1:N2
    
    threshold = thresholds(i2);
    
    tic;
    
    R1 = nonlin(U1,threshold,1);
    R2 = nonlin(U2,threshold,1);
    
    [W1,etas] = mnrfitbb(R1,Y1,W0,maxiter,conv_crit,eta0);
    ITERATIONS(i2) = length(etas);
    W0 = W1;
    
    if i2 == 1 fprintf('%d',ITERATIONS(i2)); else fprintf('-%d',ITERATIONS(i2)); end
    
    PK = mnrvals(W1,R2);
    
    FC(i2) = perfeval_fc(PK,Y2);
    
    RUNTIMES(i2) = toc/60;
end
fprintf('\n')

cd(fileparts(which(mfilename)));
clear SBANK1 SBANK2 Y1 Y2 U1 U2 R1 R2 PK

figure('units','normalized','outerposition',[0 0 1 1])
plot([])
hold on
axis([thresholds(1) thresholds(end) 0 1])
line([thresholds(1) thresholds(end)],[1/K 1/K],'linestyle',':','color','r')
xlabel('threshold [mV]');
ylabel('fraction correct');
dth = thresholds(2)-thresholds(1);
title({sprintf('N = %d, K = %d, MK = %d, J = %d, M1 = %d, M2 = %d',N,K,MK,J,M1,M2),sprintf('convcrit = %g, eta0 = %g, resolution = %g [mV]',conv_crit,eta0,dth)})
plot(thresholds,FC);

save([mfilename,'_data.mat']);
saveas(gcf,sprintf('%s_done.png',mfilename));