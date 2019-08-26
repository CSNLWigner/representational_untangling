% script for creating Fig. 4 figure supplement 1
% extrenal functions used: -
% input data:
%   MP_errorbars_data.mat
%   FR_errorbars_data.mat

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

font_size = 18;
font_ABC = 24;

bargap = 0.22;

colors = [0 0 0; 0 .3 .8; 0 .7 0; .9 0 .1; 0 .8 .9; .8 .4 1; .8 .8 .2; .7 .7 .7];
phase_color = colors(2,:);

width = 180;
height = 360;
mtop = 50;
mbottom = 65;
mright = 25;
mleft = 72;
gapx = 32;
gapy = 65;
figure_width = mleft+width+mright;
figure_height = mbottom+height+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

% panel B

load('MP_errorbars_data.mat');
FC_MP(FC_MP == .1) = NaN;
MP_MEANS = nanmean(FC_MP);
MP_ERRORS = sqrt(nanvar(FC_MP));

load('FR_errorbars_data.mat');
FR_MEANS = mean(FC_FR);
FR_ERRORS = sqrt(var(FC_FR));

axes_position = [mleft mbottom width height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

plot([])
hold on

dx = 0.25;
set(gca,'xlim',[0,8+2*dx],'xticklabels',[])
set(gca,'ylim',[0 1],'ytick',0:.1:1,'yticklabel',[])

for i1 = 1:8
    patch([i1-1+bargap/2,i1-bargap/2,i1-bargap/2,i1-1+bargap/2]+dx,[0 0 FR_MEANS(i1) FR_MEANS(i1)],colors(i1,:))
    line([i1-1+bargap/2+dx i1-bargap/2+dx],[MP_MEANS(i1) MP_MEANS(i1)],'color','w','linewidth',2)
    line([i1-0.5+dx i1-0.5+dx],[MP_MEANS(i1)-2*MP_ERRORS(i1) MP_MEANS(i1)+2*MP_ERRORS(i1)],'color','w','linewidth',2)
    line([i1-0.5+dx i1-0.5+dx],[FR_MEANS(i1) FR_MEANS(i1)+2*FR_ERRORS(i1)],'color',colors(i1,:),'linewidth',2)
end
line([1 8]+dx,[1/K 1/K],'linestyle',':','color','k')

xlabel('nuisance combination','FontName','Helvetica','fontsize',font_size)
ylab = ylabel('fraction correct','FontName','Helvetica','fontsize',font_size)
ylab_pos = get(ylab,'position');
ylab_pos(1) = ylab_pos(1)-0.5;
set(ylab,'position',ylab_pos)
set(gca,'yticklabel',{'','0.1','','','','','','','','','1'})

set(gcf,'PaperPositionMode','auto','papersize',[10 16]);
print(gcf,'Fig4_S1','-dpdf','-r0')