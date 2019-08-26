close all
clearvars

cd(fileparts(which(mfilename)));

N = 500; K = 10; MK = 10; J = 50; M2 = 5;

disk_radius = 90; gcenter_seed = 1; gphi_seed = 1;

stim_DC = 0.5; V0 = -60; V1 = 12; noise = 4.2; 

lambda0 = 3; sigma0 = 2; c0 = 0.5;

prefactor = 1; exponent = 1;

thresholds = -56:0.5:-51;

N1 = length(thresholds);

FC_OPT = zeros(1,N1);
PFC_OPT = zeros(1,N1);
RUNTIMES_OPT = zeros(1,N1);

% GBANK:

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

SBANK2 = repmat(SBANK0,M2,1);

Y2 = ceil(SBANK2(:,2)/(180/K));

MU = V0+V1*cresponse(GBANK,SBANK0); % mu_n(theta_k,phi_j) = MU((k-1)*J+j,n)

MU3D = permute(reshape(MU',N,J,K*MK),[1 3 2]); % MU3D(n,k,j) = g_n(theta_k,phi_j)

I2 = K*MK*J*M2;

U2 = repmat(MU,M2,1)+noise*randn(I2,N);

for i1 = 1:N1
    
    tic

    threshold = thresholds(i1);
    
    logrho_a = -1/noise^2/prefactor^2/2;
    logrho_b0 = -threshold/noise^2/prefactor;
    LOGRHO_B1 = MU3D/noise^2/prefactor;
    logrho_c0 = -threshold^2/noise^2/2-log(noise*prefactor*sqrt(2*pi));
    LOGRHO_C1 = -MU3D.^2/noise^2/2+threshold*MU3D/noise^2;
    P0 = .5+.5*erf((threshold-MU3D)/noise/sqrt(2));
    LOGP0 = log(P0);
    
    R2 = nonlin(U2,threshold,1);
    IS0 = (R2 == 0);
    ISNOT0 = ~IS0;
    ISZEROS = double(IS0);
    ISNOT0 = double(ISNOT0);
    
    PMK_OPT = zeros(K*MK,I2);
    for ii = 1:I2
        R2ii = R2(ii,:);
        ISZERO = IS0(ii,:);
        ISNOTZERO = ISNOT0(ii,:);
        
        ISZERO3D = repmat(ISZERO',1,K*MK,J);
        ISNOTZERO3D = repmat(ISNOTZERO',1,K*MK,J);
        R2ii3D = repmat(R2ii',1,K*MK,J);
        
        TERM1 = ISZERO3D.*LOGP0;
        TERM2 = ISNOTZERO3D.*(logrho_a*R2ii3D.^2+(logrho_b0+LOGRHO_B1).*R2ii3D+logrho_c0+LOGRHO_C1);
        
        SUMLOGii_kj = squeeze(sum(TERM1)+sum(TERM2)); % summa N cells 
        alpha0 = min(min(SUMLOGii_kj));
        
        PMKOPTii = exp(SUMLOGii_kj-alpha0);
        PMK_OPT(:,ii) = sum(PMKOPTii,2); % summa J
    end
    PMK_OPT = PMK_OPT./repmat(sum(PMK_OPT),K*MK,1); % normalization
    
    PK_OPT = zeros(K,I2);
    for kk = 1:K
        for mk = 1:MK
            PK_OPT(kk,:) = PK_OPT(kk,:)+PMK_OPT((kk-1)*MK+mk,:);
        end
    end
    
    FC_OPT(i1) = perfeval_fc(PK_OPT,Y2);
    PFC_OPT(i1) = perfeval_pfc(PK_OPT,Y2);
    
    RUNTIMES_OPT(i1) = toc;
    
    save(sprintf('%s_data.mat',mfilename),'N','K','MK','J','M2','disk_radius','gcenter_seed','gphi_seed','stim_DC','V0','V1','noise','lambda0','sigma0','c0','prefactor','exponent','thresholds','FC_OPT','PFC_OPT','RUNTIMES_OPT');
    
    fprintf('.');
end

runtime = ceil(sum(RUNTIMES_OPT(:))/60);
minute_data = mod(runtime,60);
hour_data = (runtime-minute_data)/60;
fprintf(' %dh%dm\n',hour_data,minute_data);

grey = [.8 .8 .8];
phase_color = [0 .3 .8];
light_blue = 0.33*phase_color + 0.66*grey;

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
plot([])
hold on
axis([thresholds(1) thresholds(end) 0 1])
xlabel('threshold [mV]');
ylabel('performance');
title(sprintf('N = %d, K = %d, MK = %d, J = %d, M2 = %d',N,K,MK,J,M2))
p1 = plot(thresholds,FC_OPT,'color',light_blue,'linewidth',2);
p2 = plot(thresholds,PFC_OPT,'color',light_blue,'linewidth',2,'linestyle',':');
l0 = line([thresholds(1) thresholds(end)],[1/K 1/K],'color','r','linestyle',':');
legend([p1,p2,l0],{'FC','PFC','chance level'})

saveas(gcf,sprintf('%s_done.png',mfilename));