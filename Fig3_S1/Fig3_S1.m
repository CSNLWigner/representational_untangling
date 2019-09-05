% script for creating Fig. 3 figure supplement 1
% note: running time ~ 1min
% extrenal functions used:
%   paramset
%   randisk
%   circular_gabor

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

text_size = 18;
panel_label_size = round(1.5*text_size);

width_A = 400;
width_B = 640;
gapx = 80;
gapy = 70;
subplot_height = (width_A-gapy)/2;
mtop = 55;
mbottom = 55;
mright = 82;
mleft = 28;
figure_width = mleft+width_A+gapx+width_B+mright;
figure_height = mbottom+width_A+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

% panel A

axes_position = [mleft mbottom width_A width_A];
axes('unit','pixel','position',axes_position)

N = 500;

disk_radius = 90; gcenter_seed = 1; gphi_seed = 1;

lambda0 = 3; sigma0 = 2;
  
gx0 = {'const',0};
gy0 = {'const',0};
glambda = {'const',lambda0};
gtheta = {'grid',0,180,N};
gphi = {'perm',0,360,gphi_seed};
gsigma = {'const',sigma0};
GBANK = paramset(gx0,gy0,glambda,gtheta,gphi,gsigma);
GBANK(:,1:2) = randisk(disk_radius,N,gcenter_seed);

hold on
axis square

grid_size = 1200;
F = 0;
for n = 1:N
    gabor_params = GBANK(n,:);
    G = circular_gabor([90],gabor_params,grid_size);
    F = F+G;
end

threshold = 0.02;
F(F>threshold) = threshold;
F(F<-threshold) = -threshold;

for i = 1:grid_size
    for j = 1:grid_size
        if (i-grid_size/2)^2+(j-grid_size/2)^2 > (grid_size/2)^2
            F(i,j) = threshold;
        end
    end
end

imagesc(F)
colormap(gray(256));
axis off
hold on

varphi = 0:pi/100:2*pi;
xx = grid_size/2+grid_size*cos(varphi)/2;
yy = grid_size/2+grid_size*sin(varphi)/2;
plot(xx,yy,'k','linewidth',1);

xlim([0 grid_size])
ylim([0 grid_size])

text(10,grid_size+96,'A','fontsize',panel_label_size,'FontName','Helvetica','HorizontalAlignment','center') 
    
% panel B

K = 10; MK = 10; J = 50; M1 = 10;
maxiter = 10000; conv_crit = 1e-4; eta0 = 1e-8;
stim_DC = 0.5; V0 = -60; V1 = 12; noise = 3; c0 = 0.5;

threshold_MP = -80;
threshold_FR = -57;

W_MP = zeros(N+1,K-1);
W_FR = zeros(N+1,K-1);

s0lambda = {'const',lambda0};
s0theta = {'grid',0,180,K*MK};
s0phi = {'grid',0,360,J};
s0DC = {'const',stim_DC};
s0AC = {'const',stim_DC*c0};

SBANK0 = paramset(s0lambda,s0theta,s0phi,s0DC,s0AC);
SBANK1 = repmat(SBANK0,M1,1);
Y1 = ceil(SBANK1(:,2)/(180/K));
U_MEAN = V0+V1*cresponse(GBANK,SBANK0);
U1 = repmat(U_MEAN,M1,1)+noise*randn(K*MK*J*M1,N);

threshold = threshold_MP;
tic
R1 = nonlin(U1,threshold,1);
[W1,etas] = mnrfitbb(R1,Y1,0,maxiter,conv_crit,eta0);
W_MP = W1;
fprintf('MP iterations: %d (%gsec)\n',length(etas),toc);

threshold = threshold_FR;
tic
R1 = nonlin(U1,threshold,1);
[W1,etas] = mnrfitbb(R1,Y1,0,maxiter,conv_crit,eta0);
W_FR = W1;
fprintf('FR iterations: %d (%gsec)\n',length(etas),toc);

eta_MP = W_MP(1,:);
W_MP = W_MP(2:end,:);
W_MP = [W_MP,zeros(N,1)];
W_MP = W_MP-repmat(mean(W_MP,2),1,10);
min_MP = min(min(W_MP));
max_MP = max(max(W_MP));

eta_FR = W_FR(1,:);
W_FR = W_FR(2:end,:);
W_FR = [W_FR,zeros(N,1)];
W_FR = W_FR-repmat(mean(W_FR,2),1,10);
min_FR = min(min(W_FR));
max_FR = max(max(W_FR));

axes_position = [mleft+width_A+gapx mbottom+subplot_height+gapy width_B subplot_height];
a1 = axes('unit','pixel','position',axes_position);
set(a1,'fontsize',text_size,'FontName','Helvetica')

maxcolor = [1 0 0];
midcolor = [1 1 1];
mincolor = [0 0 1];

mix = linspace(1,0,100)';
colors1 = mix*mincolor+(1-mix)*midcolor;
colors2 = mix*midcolor+(1-mix)*maxcolor;
COLORS = [colors1;colors2];

imagesc(flipud(W_MP'))

for k = 1:K
    line([1 500.5 500.5 1 1],[k-0.5 k-0.5 k+0.5 k+0.5 k-0.5],'linewidth',1,'color','k')   
end

ylab = ylabel('stimulus orientations','fontsize',text_size,'FontName','Helvetica');
ylab_pos = get(ylab,'position');
ylab_pos(1) = ylab_pos(1)-8;
set(ylab,'position',ylab_pos)

set(gca,'xtick',.5:125:500.5,'xticklabel',{sprintf('0%c',char(176)),sprintf('45%c',char(176)),sprintf('90%c',char(176)),sprintf('135%c',char(176)),sprintf('180%c',char(176))})
set(gca,'ytick',0.5:5:10.5,'yticklabel',{sprintf('180%c',char(176)),sprintf('90%c',char(176)),sprintf('0%c',char(176))})
for i = 125.5:125:500.5
    line([i i],[10.5 10.1],'linewidth',1,'color','k')   
end

title_MP = title('membrane potential-based decoding','fontsize',text_size,'FontName','Helvetica');
title_pos = get(title_MP,'position');
title_pos(2) = title_pos(2)-0.75;
set(title_MP,'position',title_pos)

cmax_MP = 0.02;
caxis([-cmax_MP cmax_MP])
colormap(a1,COLORS)
cpos = [(figure_width-0.82*mright)/figure_width (mbottom+subplot_height+gapy)/figure_height 0.015 subplot_height/figure_height];
colorbar('ytick',-0.02:0.01:0.02,'fontsize',round(0.6*text_size),'FontName','Helvetica','linewidth',1,'position',cpos)

text(-28,-1.45,'B','fontsize',panel_label_size,'FontName','Helvetica','HorizontalAlignment','center') 

% panel C

axes_position = [mleft+width_A+gapx mbottom width_B subplot_height];
a2 = axes('unit','pixel','position',axes_position);
set(a2,'fontsize',text_size,'FontName','Helvetica')

imagesc(flipud(W_FR'))

xlabel('neurons (preferred orientation)','fontsize',text_size,'FontName','Helvetica')

ylab = ylabel('stimulus orientations','fontsize',text_size,'FontName','Helvetica');
ylab_pos = get(ylab,'position');
ylab_pos(1) = ylab_pos(1)-8;
set(ylab,'position',ylab_pos)

set(gca,'xtick',.5:125:500.5,'xticklabel',{sprintf('0%c',char(176)),sprintf('45%c',char(176)),sprintf('90%c',char(176)),sprintf('135%c',char(176)),sprintf('180%c',char(176))})
set(gca,'ytick',0.5:5:10.5,'yticklabel',{sprintf('180%c',char(176)),sprintf('90%c',char(176)),sprintf('0%c',char(176))})
for i = 125.5:125:500.5
    line([i i],[10.5 10.1],'linewidth',1,'color','k')   
end

title_FR = title('firing rate-based decoding','fontsize',text_size,'FontName','Helvetica');
title_pos = get(title_FR,'position');
title_pos(2) = title_pos(2)-0.75;
set(title_FR,'position',title_pos)
for k = 1:K
    line([1 500.5 500.5 1 1],[k-0.5 k-0.5 k+0.5 k+0.5 k-0.5],'linewidth',1,'color','k')   
end
cmax_FR = 0.35;
caxis([-cmax_FR cmax_FR])
colormap(a2,COLORS)
cpos = [(figure_width-0.82*mright)/figure_width mbottom/figure_height 0.015 subplot_height/figure_height];
colorbar('ytick',-0.3:0.1:0.3,'fontsize',round(0.6*text_size),'FontName','Helvetica','linewidth',1,'position',cpos)

text(-28,-1.45,'C','fontsize',panel_label_size,'FontName','Helvetica','HorizontalAlignment','center') 

set(gcf,'PaperPositionMode','auto','papersize',[43 18]);
print(gcf,mfilename,'-dpdf','-r0')
saveas(gcf,[mfilename,'.png']);