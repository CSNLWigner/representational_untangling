% script for creating Fig. 7 figure supplement 3

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

load('kappa12_ranges_data');

font_size = 22;

color1 = [0 .3 .8];
color2 = [1 .5 0];
grey = [.84 .84 .84];

width = 310;
height = 295;
mtop = 30;
mbottom = 75;
mright = 30;
mleft = 85;
figure_width = mleft+width+mright;
figure_height = mbottom+height+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

axes_position = [mleft mbottom width height];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

plot([])
hold on

xrange = [0 .6];
yrange = [-.6 1];

set(gca,'FontName','Helvetica','fontsize',font_size)
set(gca,'xlim',xrange,'ylim',yrange,'fontsize',font_size,'linewidth',1.5)
set(gca,'xtick',.2:.2:.6,'ytick',-.6:.2:1)

V0 = -60; V1 = 12;

sigma_values = 1.5:0.5:7;

y1 = OPTTH(1,:);
y1_dn = THMIN90(1,:);
y1_up = THMAX90(1,:);

y2 = OPTTH(2,:);
y2_dn = THMIN90(2,:);
y2_up = THMAX90(2,:);

f1 = fill([sigma_values,fliplr(sigma_values)]/V1,[y1_up-V0,fliplr(y1_dn)-V0]/V1,color1,'edgecolor','none')
set(f1,'facealpha',.2)
f2 = fill([sigma_values,fliplr(sigma_values)]/V1,[y2_up-V0,fliplr(y2_dn)-V0]/V1,color2,'edgecolor','none')
set(f2,'facealpha',.2)

plot(sigma_values/V1,(y1-V0)/V1,'color',color1,'linewidth',2)
plot(sigma_values/V1,(y2-V0)/V1,'color',color2,'linewidth',2)

set(gca,'layer','top')

xlabel('noise-to-signal ratio','FontName','Helvetica','fontsize',font_size)
ylabel('normalised threshold','FontName','Helvetica','fontsize',font_size)

set(gcf,'PaperPositionMode','auto','papersize',[15.5 14.5]);
print(gcf,mfilename,'-dpdf','-r0')
saveas(gcf,[mfilename,'.png']);