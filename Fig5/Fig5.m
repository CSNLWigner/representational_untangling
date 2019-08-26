% script for creating Fig. 5
% input data:
%   script_4d_lognorm_data.mat
%   script_4d_beta_data.mat

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

font_size = 18;
bigfont_size = 24;

light_grey = [.9 .9 .9];

width = 320;
height1 = 300;
height2 = 120;
mtop = 35;
mbottom = 45;
mright = 20;
mleft = 70;
gapx = 50;
gapy = 85;
barwidth = 15;
smallbargap = 2;

figure_width = mleft+2*width+gapx+mright;
figure_height = mbottom+height2+gapy+height1+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

% panel A

axes_position = [mleft mbottom+height2+gapy width height1];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

text(-0.18,1.05,'A','fontname','Helvetica','fontsize',bigfont_size,'units','normalized')
plot([])
hold on

xmin = 0.01; xmax = 20; dx = 0.01;
pmu_values = [0.5 0.95 1.4 1.85];
psigma_values = [0.55 0.55 0.55 0.55];
colors1 = [0 .8 0; 0 .7 0; 0 .6 0; 0 .5 0];

m = 2;
pmu = pmu_values(m);
psigma = psigma_values(m);
lognorm = @(x) exp(-(log(x)-pmu).^2/psigma^2/2)./x/psigma/sqrt(2*pi);
xx = xmin:dx:xmax;
yy = lognorm(xx)/integral(lognorm,xmin,xmax);
area(xx,yy,'FaceColor',light_grey)
set(gca,'layer','top','linewidth',1) 

M1 = length(pmu_values);
for m = 1:M1
    pmu = pmu_values(m);
    psigma = psigma_values(m);
    lognorm = @(x) exp(-(log(x)-pmu).^2/psigma^2/2)./x/psigma/sqrt(2*pi);
    xx = xmin:dx:xmax;
    yy = lognorm(xx)/integral(lognorm,xmin,xmax);
    plot(xx,yy,'color',colors1(m,:),'linewidth',2)
end

set(gca,'xlim',[0 xmax],'ylim',[0 0.55])
xlabel('spatial period [deg]','fontname','Helvetica','fontsize',font_size)
ylabel('probability density','fontname','Helvetica','fontsize',font_size)
ylabh = get(gca,'ylabel');
set(ylabh,'position',get(ylabh,'position') - [0.45 0 0])

patch([5.5 7.5 7.5 5.5],[0.47 0.47 0.53 0.53],light_grey)
text(8.5,0.5,'reference distribution','fontname','Helvetica','fontsize',font_size)

% panel B

axes_position = [mleft+width+gapx mbottom+height2+gapy width height1];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

text(-0.1,1.05,'C','FontName','Helvetica','fontsize',bigfont_size,'units','normalized')
plot([])
hold on

c_resolution = 0.01;
palpha_values = [1.7 2.4 3.6 3.6 3.6];
pbeta_values = [3.6 3.6 3.6 2.4 1.7];
colors2 = [1 0 0; .9 0 .1; .8 0 .2; .7 0 .3; .6 0 .4];

m = 2;
palpha = palpha_values(m);
pbeta = pbeta_values(m);
beta = @(x) (x.^(palpha-1)).*(1-x).^(pbeta-1);
xx = 0:c_resolution:1;
yy = beta(xx)/integral(beta,0,1);
area(xx,yy,'FaceColor',[.9 .9 .9])
set(gca,'layer','top','linewidth',1) 

M2 = length(palpha_values);
for m = 1:M2
    palpha = palpha_values(m);
    pbeta = pbeta_values(m);
    beta = @(x) (x.^(palpha-1)).*(1-x).^(pbeta-1);
    xx = 0:c_resolution:1;
    yy = beta(xx)/integral(beta,0,1);
    plot(xx,yy,'color',colors2(m,:),'linewidth',2)
end

set(gca,'xlim',[0 1],'fontsize',font_size)
set(gca,'ylim',[0 2.25],'fontsize',font_size)
xlabel('contrast','FontName','Helvetica','fontsize',font_size)

% panel C

axes_position = [mleft mbottom width height2];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

text(-0.185,1.25,'B','fontsize',bigfont_size,'units','normalized')
plot([])
hold on

set(gca,'xlim',[0 M1+1],'xtick',[])
set(gca,'ylim',[0 .5],'ytick',[0 .2 .4])
xlabel('test distribution','FontName','Helvetica','fontsize',font_size)
ylabel('fraction correct','FontName','Helvetica','fontsize',font_size)
ylabh = get(gca,'ylabel');
set(ylabh,'position',get(ylabh,'position') - [0.1 0 0])
 
load('script_4d_lognorm_data.mat','OPT_FC','GREY_FC')
 
dm = 0.3;
order = [2 1 3 4];
for m = 1:M1
       line([m m],[0 GREY_FC(order(m))],'color',colors1(m,:),'linewidth',25)
       if m == 2
           line([m m],[0 GREY_FC(order(m))],'color',[.9 .9 .9],'linewidth',20)
       end
       line([m-dm m+dm],[OPT_FC(order(m)) GREY_FC(order(m))],'color',colors1(m,:),'linewidth',2)
end
set(gca,'layer','top') 

% panel D

axes_position = [mleft+width+gapx mbottom width height2];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

text(-0.1,1.25,'D','fontsize',bigfont_size,'units','normalized')
plot([])
hold on

set(gca,'xlim',[0 M2+1],'xtick',[])
set(gca,'ylim',[0 .5],'ytick',[0 .2 .4])
xlabel('test distribution','FontName','Helvetica','fontsize',font_size)

load('script_4d_beta_data.mat','OPT_FC','GREY_FC')

dm = 0.35;
for m = 1:M2
      line([m m],[0 GREY_FC(m)],'color',colors2(m,:),'linewidth',25)
      if m == 2
          line([m m],[0 GREY_FC(m)],'color',[.9 .9 .9],'linewidth',20)
      end
      line([m-dm m+dm],[OPT_FC(m) OPT_FC(m)],'color',colors2(m,:),'linewidth',2)
end
set(gca,'layer','top') 

set(gcf,'PaperPositionMode','auto','papersize',[28 21]);
print(gcf,'Fig5','-dpdf','-r0')
saveas(gcf,sprintf('%s.png',mfilename));