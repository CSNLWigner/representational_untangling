% script for creating Fig. 3 figure supplement 2
% extrenal functions used: -
% input data:
%   script_2d_3deg_data.mat (3deg curves)
%   script_2d_90deg_data.mat (90deg curves)

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

load('script_2d_3deg_data.mat','thresholds','FC_LIN')
FC = FC_LIN;
FC_min = min(FC);
FC_max = max(FC);

load('script_2d_90deg_data.mat','thresholds','FC_LIN')
FC90 = FC_LIN;
FC90_min = min(FC90);
FC90_max = max(FC90);

text_size = 16;
ABC_size = 24;
ABC_yshift = 0.1;

lwidth = 3;
ax_width = 1;
xmin = -72.5;
xmax = -41;
ymax = 1;

width_A = 400;
height = 320;
mtop = 25;
mbottom = 60;
mright = 27;
mleft = 67;
gapx = 60;

figure_width = mleft+width_A+mright;
figure_height = mbottom+height+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

% panel A

axes_position = [mleft mbottom width_A height];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',text_size,'linewidth',ax_width)

hold on

f1 = fill([thresholds,fliplr(thresholds)],[FC_max,fliplr(FC_min)],'r','edgecolor','none')
set(f1,'facealpha',.4)
f2 = fill([thresholds,fliplr(thresholds)],[FC90_max,fliplr(FC90_min)],'b','edgecolor','none')
set(f2,'facealpha',.4)

K = 10;
line([thresholds(1) thresholds(end)],[1/K 1/K],'linestyle',':','linewidth',ax_width,'color','k')

set(gca,'xlim',[xmin,xmax],'xtick',-70:5:-45,'xticklabel',{'-70','-65','-60','-55','-50','-45'})
set(gca,'ylim',[0 ymax],'ytick',0:.1:ymax)

xlabel('FRNL threshold [mV]','FontName','Helvetica','fontsize',text_size)
ylabel('fraction correct','FontName','Helvetica','fontsize',text_size)

set(gcf,'PaperPositionMode','auto','papersize',[17 14]);
print(gcf,mfilename,'-dpdf','-r0')
saveas(gcf,[mfilename,'.png']);