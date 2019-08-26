% script for creating Fig. 3
% extrenal functions used:
%   cresponse
%   circular_gabor
%   paramset
%   maxfit
% input data:
%   script_2d_data.mat (performance curve of the linear decoder: phase nuisance)
%   script_1d_data.mat (performance curve of the linear decoder: no nuisance)
%   script_2d_Bayes_data.mat (performance curve of the optilmal decoder: phase nuisance)
%   script_1d_Bayes_data.mat (performance curve of the optilmal decoder: no nuisance)
%   script_2d_sparseness_data.mat (sparseness curve)
%   histogram_data.mat (MP histogram)
%   script_1d_pixel_DoG_Gabor_data.mat (decoding from pixels and DoG filters: no nuisance)
%   script_2d_pixel_DoG_Gabor_data.mat (decoding from pixels and DoG filters: phase nuisance)

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

load('script_2d_data.mat','thresholds','FC_LIN','K');
LINPERF_2D = FC_LIN;

load('script_1d_data.mat','FC_LIN');
LINPERF_1D = FC_LIN;

load('script_2d_Bayes_data.mat','FC_OPT');
OPTPERF_2D = FC_OPT;

load('script_1d_Bayes_data.mat','FC_OPT');
OPTPERF_1D = FC_OPT;

load('script_2d_sparseness_data.mat','SPARSENESS');
SPARSENESS(isnan(SPARSENESS)) = 1;

% hold on
% plot(0:0.01:1,0:0.01:1,'k:')
% plot(SPARSENESS.*(OPTPERF_2D-1/K),(LINPERF_2D-1/K),'linewidth',2)
% set(gca,'xlim',[0 1],'ylim',[0 1])
% axis square
% xlabel('SPARSENESS*(OPTPERF2D-1/K)')
% ylabel('LINPERF2D-1/K')
% corr = corrcoef(SPARSENESS.*(OPTPERF_2D-1/K),(LINPERF_2D-1/K));
% title(sprintf('non-diagonal corrcoef = %g',corr(1,2)))
% saveas(gcf,'corrcoef_sparseness_estimation.png');

font_size = 16;
ABC_size = 24;
leg_size = 16;
curve_width = 3;
ax_width = 1;
ABC_yshift = 0.1;
MPFR_yshift = -0.045;

xmin = -74;
xmax = -41;

% predefined colors:
phase_color = [0 .3 .8];
white = [1 1 1];
light_blue = 0.33*phase_color + 0.66*white;
grey = .75*white;
sparseness_color = [255 135 137]/255;
hist_color = [0.66 0.8 0.66];
light_grey = .9*white;

% figure settings:
width_A = 420;
height = 320;
mtop = 57;
mbottom = 60;
mright = 25;
mleft = 68;
gapx = 128;
barwidth = 25;
bargap0 = 16;
bargap1 = 6;
bargap2 = 15;
bargap3 = 14;
width_B = bargap0+4*barwidth+2*bargap1+bargap2+bargap3;
gapx2 = 40;

figure_width = mleft+width_A+gapx+3*width_B+2*gapx2+mright;
figure_height = mbottom+height+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

% panel A

axes_position = [mleft mbottom width_A height];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',ax_width)

plot([])
hold on

% generating histogram data:
% M = 50;
% load('script_2d_sparseness_data.mat','N','K','MK','J','lambda0','c0','stim_DC','V0','V1','noise','GBANK');
% s0lambda = {'const',lambda0};
% s0theta = {'grid',0,180,K*MK};
% s0phi = {'grid',0,360,J};
% s0DC = {'const',stim_DC};
% s0AC = {'const',stim_DC*c0};
% SBANK0 = paramset(s0lambda,s0theta,s0phi,s0DC,s0AC);
% SBANK = repmat(SBANK0,M,1);
% U = V0+V1*cresponse(GBANK,SBANK)+noise*randn(K*MK*J*M,N);
% [hist_counts,hist_centers] = hist(U(:),100);
% save('histogram_data.mat','hist_counts','hist_centers')
% clear SBANK0 SBANK U

load('histogram_data.mat','hist_counts','hist_centers')
area(hist_centers,0.45*hist_counts/max(hist_counts),'faceColor',hist_color,'edgecolor','w')    
set(gca,'layer','top')

plot(thresholds,SPARSENESS,'linewidth',curve_width,'color',sparseness_color)

set(gca,'xlim',[xmin xmax],'xtick',-70:10:-45)
set(gca,'ylim',[0 1],'ytick',0:.1:1,'yticklabel',{'0','','0.2','','0.4','','0.6','','0.8','','1'})

xlabel('FRNL threshold [mV]    ','FontName','Helvetica','fontsize',font_size)
ylab = ylabel('fraction correct /\color[rgb]{1 0.5294 0.5373}sparseness','FontName','Helvetica','fontsize',font_size);
ylab_pos = get(ylab,'position');
ylab_pos(1) = ylab_pos(1)-0.2;
set(ylab,'position',ylab_pos)

plot(thresholds,OPTPERF_1D,'linewidth',curve_width,'color',grey);
plot(thresholds,OPTPERF_2D,'linewidth',curve_width,'color',light_blue);
plot(thresholds,LINPERF_1D,'linewidth',curve_width,'color','k');
plot(thresholds,LINPERF_2D,'linewidth',curve_width,'color',phase_color);

line([thresholds(1) thresholds(end)],[1/K 1/K],'linestyle',':','linewidth',ax_width,'color','k')

leg_x0 = -43.5;
leg_width = 2;
leg_gapx = 3;
leg_y0 = 0.58;
leg_gapy = 0.1;
leg_titlegap = 0.09;
leg_linespace = 0.065;
leg_neg = 1.4;

text(leg_x0+0.5*leg_width,leg_y0+leg_gapy+leg_titlegap+leg_linespace,'fixed','fontsize',leg_size,'HorizontalAlignment','center','fontname','Helvetica')
text(leg_x0+1.5*leg_width+leg_gapx,leg_y0+leg_gapy+leg_titlegap+leg_linespace,'variable','fontsize',leg_size,'HorizontalAlignment','center','fontname','Helvetica')
text(leg_x0+0.5*leg_width,leg_y0+leg_gapy+leg_titlegap,'phase','fontsize',leg_size,'HorizontalAlignment','center','fontname','Helvetica')
text(leg_x0+1.5*leg_width+leg_gapx,leg_y0+leg_gapy+leg_titlegap,'phase','fontsize',leg_size,'HorizontalAlignment','center','fontname','Helvetica')

text(leg_x0-leg_neg,leg_y0+leg_gapy,'optimal','fontsize',leg_size,'HorizontalAlignment','right','fontname','Helvetica')
text(leg_x0-leg_neg,leg_y0,'linear','fontsize',leg_size,'HorizontalAlignment','right','fontname','Helvetica')

x1 = 0.3802;
x2 = 0.4328;
xw = 0.024;
y1 = 0.560;
y2 = 0.635;
annotation('line',[x1 x1+xw],[y2 y2],'linewidth',curve_width,'color',grey)
annotation('line',[x1 x1+xw],[y1 y1],'linewidth',curve_width,'color','k')
annotation('line',[x2 x2+xw],[y2 y2],'linewidth',curve_width,'color',light_blue)
annotation('line',[x2 x2+xw],[y1 y1],'linewidth',curve_width,'color',phase_color)

barheight_1d_MP = LINPERF_1D(1);
[~,maxi] = max(LINPERF_1D);
[maxx,maxy] = maxfit(thresholds,LINPERF_1D,3);
barheight_1d_FR = maxy;
Bayesbound_1d_MP = OPTPERF_1D(1);
Bayesbound_1d_FR = OPTPERF_1D(maxi);
scatter(maxx,maxy,100,'k','d','filled')

barheight_2d_MP = LINPERF_2D(1);
[~,maxi] = max(LINPERF_2D);
[maxx,maxy] = maxfit(thresholds,LINPERF_2D,3);
barheight_2d_FR = maxy;
Bayesbound_2d_MP = OPTPERF_2D(1);
Bayesbound_2d_FR = OPTPERF_2D(maxi);
scatter(maxx,maxy,100,phase_color,'d','filled')

text(-75.8,1+ABC_yshift,'A','fontsize',ABC_size,'HorizontalAlignment','center','fontweight','bold')

gwidth = 0.25*width_B;
axes_position = [mleft+width_A+gapx-2.5*gwidth mbottom+height-.45*gwidth gwidth gwidth];
subplot_axes = axes('unit','pixel','position',axes_position);
hold on
axis square
axis off

gabor_params = [0,0,3,45,0,2];

% Gabor patch for illustration
NN = 1000;
G1 = circular_gabor([5],gabor_params,NN);
maxg = max(max(G1));
for i = 1:NN
    for j = 1:NN
        if (i-NN/2)^2+(j-NN/2)^2 > (NN/2)^2
            G1(i,j) = maxg;
        end
    end 
end  
imagesc(G1)
colormap(0.25+0.75*gray(256));
cx0 = NN/2;
cy0 = NN/2;
rr = NN/2-1;
tt = linspace(0,2*pi);
xx = cx0+rr*cos(tt);
yy = cy0+rr*sin(tt);
hold on
plot(xx,yy,'color','k','linewidth',1)

% inset

axes_position = [mleft+0.085*width_A mbottom+0.50*height 0.20*width_A 0.22*height];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',10,'linewidth',1)

hold on

plot(0:0.01:1,0:0.01:1,'k')
scatter(LINPERF_2D,1/K+SPARSENESS.*(OPTPERF_2D-1/K),5,phase_color,'filled')
set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',[0 1],'ytick',[0 1],'fontname','Helvetica','fontsize',10)
axis square
xlabel({'actual performance'},'fontname','Helvetica','fontsize',10)
ylabel({'predicted performance'},'fontname','Helvetica','fontsize',10)
title({'linear decoder',''})

% panel B

load('script_1d_pixel_DoG_Gabor_data.mat','FC_Gabor')
barheight_1d_MP = FC_Gabor(1);
barheight_1d_FR = max(FC_Gabor);
load('script_2d_pixel_DoG_Gabor_data.mat','FC_Gabor')
barheight_2d_MP = FC_Gabor(1);
barheight_2d_FR = max(FC_Gabor);

axes_position = [mleft+width_A+gapx mbottom width_B height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',ax_width)

plot([])
hold on

set(gca,'xlim',[0,width_B],'xtick',[])
set(gca,'ylim',[0 1],'ytick',0:.1:1,'yticklabel',{'0','','','','','','','','','','1'})
ylab = ylabel('fraction correct','FontName','Helvetica','fontsize',font_size);
ylab_pos = get(ylab,'position');
ylab_pos(1) = ylab_pos(1)+8;
set(ylab,'position',ylab_pos)

area([bargap0 bargap0 bargap0+barwidth bargap0+barwidth],[0 barheight_1d_MP barheight_1d_MP 0],'FaceColor',light_grey)
area([bargap0+barwidth+bargap1 bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1 bargap0+2*barwidth+bargap1],[0 barheight_1d_FR barheight_1d_FR 0],'FaceColor',light_grey)
plot([bargap0 bargap0 bargap0+barwidth bargap0+barwidth],[0 barheight_1d_MP barheight_1d_MP 0],'k','LineWidth',curve_width)
plot([bargap0+barwidth+bargap1 bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1 bargap0+2*barwidth+bargap1],[0 barheight_1d_FR barheight_1d_FR 0],'k','LineWidth',curve_width)
scatter(bargap0+1.5*barwidth+bargap1,barheight_1d_FR,100,'k','d','filled')
line([bargap0 bargap0+barwidth],[Bayesbound_1d_MP Bayesbound_1d_MP],'LineWidth',curve_width,'color',grey)
line([bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1],[Bayesbound_1d_FR Bayesbound_1d_FR],'LineWidth',curve_width,'color',grey)
 
text(bargap0+barwidth/2,MPFR_yshift,'MP','FontName','Helvetica','fontsize',font_size,'HorizontalAlignment','center')
text(bargap0+bargap1+1.5*barwidth,MPFR_yshift,'FR','FontName','Helvetica','fontsize',font_size,'HorizontalAlignment','center')

bar_shift = +bargap2+bargap1+2*barwidth;
area([bargap0 bargap0 bargap0+barwidth bargap0+barwidth]+bargap2+bargap1+2*barwidth,[0 barheight_2d_MP barheight_2d_MP 0],'FaceColor',light_grey)
area([bargap0+barwidth+bargap1 bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1 bargap0+2*barwidth+bargap1]+bar_shift,[0 barheight_2d_FR barheight_2d_FR 0],'FaceColor',light_grey)
plot([bargap0 bargap0 bargap0+barwidth bargap0+barwidth]+bargap2+bargap1+2*barwidth,[0 barheight_2d_MP barheight_2d_MP 0],'color',phase_color,'LineWidth',curve_width)
plot([bargap0+barwidth+bargap1 bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1 bargap0+2*barwidth+bargap1]+bar_shift,[0 barheight_2d_FR barheight_2d_FR 0],'color',phase_color,'LineWidth',curve_width)
scatter(bargap0+1.5*barwidth+bargap1+bar_shift,barheight_2d_FR,100,phase_color,'d','filled')
line([bargap0 bargap0+barwidth]+bargap2+bargap1+2*barwidth,[Bayesbound_2d_MP Bayesbound_2d_MP],'LineWidth',curve_width,'color',light_blue)
line([bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1]+bargap2+bargap1+2*barwidth,[Bayesbound_2d_FR Bayesbound_2d_FR],'LineWidth',curve_width,'color',light_blue)
 
text(bargap0+barwidth/2+bargap2+bargap1+2*barwidth,MPFR_yshift,'MP','FontName','Helvetica','fontsize',font_size,'HorizontalAlignment','center')
text(bargap0+bargap1+1.5*barwidth+bargap2+bargap1+2*barwidth,MPFR_yshift,'FR','FontName','Helvetica','fontsize',font_size,'HorizontalAlignment','center')
 
line([0 width_B],[1/K 1/K],'linestyle',':','linewidth',ax_width,'color','k')
 
text(-20,1+ABC_yshift,'B','fontsize',ABC_size,'HorizontalAlignment','center','fontweight','bold')

% panel C

load('script_1d_pixel_DoG_Gabor_data.mat','FC_Gabor')
[~,maxindex] = max(FC_Gabor);
load('script_1d_pixel_DoG_Gabor_data.mat','FC_pixel')
barheight_1d_MP = FC_pixel(1);
barheight_1d_FR = FC_pixel(maxindex);
barheight_2d_MP = 1/K;
barheight_2d_FR = 1/K;

axes_position = [mleft+width_A+gapx+width_B+gapx2 mbottom width_B height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',ax_width)

plot([])
hold on

set(gca,'xlim',[0,width_B],'xtick',[])
set(gca,'ylim',[0 1],'ytick',0:.1:1,'yticklabel',{'0','','','','','','','','','','1'})

area([bargap0 bargap0 bargap0+barwidth bargap0+barwidth],[0 barheight_1d_MP barheight_1d_MP 0],'FaceColor',light_grey)
area([bargap0+barwidth+bargap1 bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1 bargap0+2*barwidth+bargap1],[0 barheight_1d_FR barheight_1d_FR 0],'FaceColor',light_grey)
plot([bargap0 bargap0 bargap0+barwidth bargap0+barwidth],[0 barheight_1d_MP barheight_1d_MP 0],'k','LineWidth',curve_width)
plot([bargap0+barwidth+bargap1 bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1 bargap0+2*barwidth+bargap1],[0 barheight_1d_FR barheight_1d_FR 0],'k','LineWidth',curve_width)

text(bargap0+barwidth/2,MPFR_yshift,'MP','FontName','Helvetica','fontsize',font_size,'HorizontalAlignment','center')
text(bargap0+bargap1+1.5*barwidth,MPFR_yshift,'FR','FontName','Helvetica','fontsize',font_size,'HorizontalAlignment','center')

bar_shift = +bargap2+bargap1+2*barwidth;
area([bargap0 bargap0 bargap0+barwidth bargap0+barwidth]+bargap2+bargap1+2*barwidth,[0 barheight_2d_MP barheight_2d_MP 0],'FaceColor',light_grey)
area([bargap0+barwidth+bargap1 bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1 bargap0+2*barwidth+bargap1]+bar_shift,[0 barheight_2d_FR barheight_2d_FR 0],'FaceColor',light_grey)
plot([bargap0 bargap0 bargap0+barwidth bargap0+barwidth]+bargap2+bargap1+2*barwidth,[0 barheight_2d_MP barheight_2d_MP 0],'color',phase_color,'LineWidth',curve_width)
plot([bargap0+barwidth+bargap1 bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1 bargap0+2*barwidth+bargap1]+bar_shift,[0 barheight_2d_FR barheight_2d_FR 0],'color',phase_color,'LineWidth',curve_width)

text(bargap0+barwidth/2+bargap2+bargap1+2*barwidth,MPFR_yshift,'MP','FontName','Helvetica','fontsize',font_size,'HorizontalAlignment','center')
text(bargap0+bargap1+1.5*barwidth+bargap2+bargap1+2*barwidth,MPFR_yshift,'FR','FontName','Helvetica','fontsize',font_size,'HorizontalAlignment','center')
 
line([0 width_B],[1/K 1/K],'linestyle',':','linewidth',ax_width,'color','k')
 
text(-17,1+ABC_yshift,'C','fontsize',ABC_size,'HorizontalAlignment','center','fontweight','bold')

gwidth = 0.25*width_B;
axes_position = [mleft+width_A+gapx+width_B+gapx2+1.7*gwidth mbottom+height-.45*gwidth gwidth gwidth];
subplot_axes = axes('unit','pixel','position',axes_position);
hold on
axis square
axis off

% filter illustration
NN = 1000;
G1 = zeros(NN,NN);
G1(500,500) = -1;
maxg = max(max(G1));
for i = 1:NN
    for j = 1:NN
        if (i-NN/2)^2+(j-NN/2)^2 > (NN/2)^2
            G1(i,j) = 1;
        end
    end 
end  
imagesc(G1)
colormap(0.25+0.75*gray(256));
cx0 = NN/2;
cy0 = NN/2;
rr = NN/2-1;
tt = linspace(0,2*pi);
xx = cx0+rr*cos(tt);
yy = cy0+rr*sin(tt);
hold on
plot(xx,yy,'color','k','linewidth',1)

scatter(cx0,cy0,20,'w','filled')

% panel D

load('script_1d_pixel_DoG_Gabor_data.mat','FC_Gabor')
[~,maxindex] = max(FC_Gabor);
load('script_1d_pixel_DoG_Gabor_data.mat','FC_DoG')
barheight_1d_MP = FC_DoG(1);
barheight_1d_FR = FC_DoG(maxindex);
barheight_2d_MP = 1/K;
barheight_2d_FR = 1/K;

axes_position = [mleft+width_A+gapx+2*width_B+2*gapx2 mbottom width_B height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',ax_width)

plot([])
hold on

set(gca,'xlim',[0,width_B],'xtick',[])
set(gca,'ylim',[0 1],'ytick',0:.1:1,'yticklabel',{'0','','','','','','','','','','1'})

area([bargap0 bargap0 bargap0+barwidth bargap0+barwidth],[0 barheight_1d_MP barheight_1d_MP 0],'FaceColor',light_grey)
area([bargap0+barwidth+bargap1 bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1 bargap0+2*barwidth+bargap1],[0 barheight_1d_FR barheight_1d_FR 0],'FaceColor',light_grey)
plot([bargap0 bargap0 bargap0+barwidth bargap0+barwidth],[0 barheight_1d_MP barheight_1d_MP 0],'k','LineWidth',curve_width)
plot([bargap0+barwidth+bargap1 bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1 bargap0+2*barwidth+bargap1],[0 barheight_1d_FR barheight_1d_FR 0],'k','LineWidth',curve_width)

text(bargap0+barwidth/2,MPFR_yshift,'MP','FontName','Helvetica','fontsize',font_size,'HorizontalAlignment','center')
text(bargap0+bargap1+1.5*barwidth,MPFR_yshift,'FR','FontName','Helvetica','fontsize',font_size,'HorizontalAlignment','center')

bar_shift = +bargap2+bargap1+2*barwidth;
area([bargap0 bargap0 bargap0+barwidth bargap0+barwidth]+bargap2+bargap1+2*barwidth,[0 barheight_2d_MP barheight_2d_MP 0],'FaceColor',light_grey)
area([bargap0+barwidth+bargap1 bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1 bargap0+2*barwidth+bargap1]+bar_shift,[0 barheight_2d_FR barheight_2d_FR 0],'FaceColor',light_grey)
plot([bargap0 bargap0 bargap0+barwidth bargap0+barwidth]+bargap2+bargap1+2*barwidth,[0 barheight_2d_MP barheight_2d_MP 0],'color',phase_color,'LineWidth',curve_width)
plot([bargap0+barwidth+bargap1 bargap0+barwidth+bargap1 bargap0+2*barwidth+bargap1 bargap0+2*barwidth+bargap1]+bar_shift,[0 barheight_2d_FR barheight_2d_FR 0],'color',phase_color,'LineWidth',curve_width)

text(bargap0+barwidth/2+bargap2+bargap1+2*barwidth,MPFR_yshift,'MP','FontName','Helvetica','fontsize',font_size,'HorizontalAlignment','center')
text(bargap0+bargap1+1.5*barwidth+bargap2+bargap1+2*barwidth,MPFR_yshift,'FR','FontName','Helvetica','fontsize',font_size,'HorizontalAlignment','center')
 
line([0 width_B],[1/K 1/K],'linestyle',':','linewidth',ax_width,'color','k')
 
text(-17,1+ABC_yshift,'D','fontsize',ABC_size,'HorizontalAlignment','center','fontweight','bold')

gwidth = 0.25*width_B;
axes_position = [mleft+width_A+gapx+2*width_B+2*gapx2+1.7*gwidth mbottom+height-.45*gwidth gwidth gwidth];
subplot_axes = axes('unit','pixel','position',axes_position);
hold on
axis square
axis off

gabor_params1 = [0,0,1000,0,0,.75];
gabor_params2 = [0,0,1000,0,0,2];

% DoG filter patch for illustration
NN = 1000;
G1 = .4+.2*circular_gabor([5],gabor_params1,NN)-circular_gabor([5],gabor_params2,NN);
maxg = max(max(G1));
for i = 1:NN
    for j = 1:NN
        if (i-NN/2)^2+(j-NN/2)^2 > (NN/2)^2
            G1(i,j) = maxg;
        end
    end 
end  
imagesc(G1)
colormap(0.25+0.75*gray(256));
cx0 = NN/2;
cy0 = NN/2;
rr = NN/2-1;
tt = linspace(0,2*pi);
xx = cx0+rr*cos(tt);
yy = cy0+rr*sin(tt);
hold on
plot(xx,yy,'color','k','linewidth',1)

set(gcf,'PaperPositionMode','auto','papersize',[42 15]);
print(gcf,'Fig3','-dpdf','-r0')
saveas(gcf,sprintf('%s.png',mfilename));