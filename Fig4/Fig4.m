% script for creating Fig. 4
% extrenal functions used: maxfit
% input data: script_1d2d3d4d_data.mat

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

text_size = 18;
panel_label_size = round(1.5*text_size);
lwidth = 3;
bargap = 0.22;
panel_label_yshift = 0.08;

colors = [0 0 0; 0 .3 .8; 0 .7 0; .9 0 .1; 0 .8 .9; .8 .4 1; .8 .8 .2; .7 .7 .7];

xmin = -72.5;
xmax = -41;

width_A = 560;
width_C = 180;
height_A = 120;
height_B = 360;
mtop = 52;
mbottom = 62;
mright = 25;
mleft = 70;
gapx = 36;
gapy = 95;
figure_width = mleft+width_A+gapx+width_C+mright;
figure_height = mbottom+height_B+gapy+height_A+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')
subgapx = 38;
subwidth = (width_A+gapx+width_C-2*subgapx)/3;

load('script_1d2d3d4d_data.mat');

nuisance(1).vrange = [0 360];
nuisance(1).star = 180;

nuisance(2).vrange = [0 8];
nuisance(2).star = 3;

nuisance(3).vrange = [0 1];
nuisance(3).star = 0.5;

% panel A

for i1 = 1:3
    axes_position = [mleft+(i1-1)*(subwidth+subgapx) mbottom+height_B+gapy subwidth height_A];
    axes('unit','pixel','position',axes_position)
    set(gca,'FontName','Helvetica','fontsize',text_size,'linewidth',1)
    
    hold on
    box off
    
    set(gca,'ylim',[0 2],'ytick',[],'yticklabels',[],'fontsize',text_size)
    set(gca,'xlim',nuisance(i1).vrange,'xtick',nuisance(i1).range,'fontsize',text_size)
    
    xx = linspace(nuisance(i1).range(1),nuisance(i1).range(2),50);
    switch i1
        case 1
            yy = 0*xx+0.75;
            xlab = xlabel('phase','FontName','Helvetica','fontsize',text_size)
            set(gca,'xticklabel',{sprintf('0%c',char(176)),sprintf('360%c',char(176))})
        case 2
            yy = 5*lognorm(xx)/integral(lognorm,0.5,7.5);
            xlab = xlabel('spatial period','FontName','Helvetica','fontsize',text_size)
            set(gca,'xticklabel',{sprintf('0.5%c',char(176)),sprintf('7.5%c',char(176))})
        case 3
            yy = 0.6*beta(xx)/integral(beta,0,1);
            xlab = xlabel('contrast','FontName','Helvetica','fontsize',text_size)
    end
    
    xlab_pos = get(xlab,'position');
    xlab_pos(2) = xlab_pos(2)+0.25;
    set(xlab,'position',xlab_pos)
    
    plot([])
    area(xx,yy,'FaceColor',[.9 .9 .9],'EdgeColor','none')
    set(gca,'layer','top')
    plot(xx,yy,'color',colors(i1+1,:),'linewidth',lwidth)
    plot([nuisance(i1).star nuisance(i1).star],[0 0.1],'k','linewidth',2)
   
    if i1 == 1
        text(-65,0.43,'probability','FontName','Helvetica','fontsize',text_size,'rotation',90)
        text(-35,0.63,'density','FontName','Helvetica','fontsize',text_size,'rotation',90)
        text(-70,2.5,'A','fontsize',panel_label_size,'FontName','Helvetica','HorizontalAlignment','center')
    end     
end

% panel B

axes_position = [mleft mbottom width_A height_B];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',text_size,'linewidth',1)

text(-75.1,1+panel_label_yshift,'B','fontsize',panel_label_size,'FontName','Helvetica','HorizontalAlignment','center')
hold on

set(gca,'xlim',[xmin,xmax],'xtick',-70:5:-45,'xticklabel',{'-70','-65','-60','-55','-50','-45'})
set(gca,'ylim',[0 1],'ytick',0:.1:1)
xlabel('FRNL threshold [mV]','FontName','Helvetica','fontsize',text_size)
ylab = ylabel('fraction correct','FontName','Helvetica','fontsize',text_size)

ylab_pos = get(ylab,'position');
ylab_pos(1) = ylab_pos(1)-0.15;
set(ylab,'position',ylab_pos)

plots = [];

leg_x0 = -47.5;
leg_y0 = 0.64;
leg_width = 1.75;
leg_dx = 1.4;
leg_dy = 0.05;
leg_xgap = 1.4;
leg_ygap = 0.035;
leg_yshift = 0.005;
dot_yshift = 0.01;
leg_xshift = -0.11;

ch0 = '.';
ch1 = 'v';
legends = {ch0 ch0 ch0; ch1 ch0 ch0; ch0 ch1 ch0; ch0 ch0 ch1; ch1 ch1 ch0; ch1 ch0 ch1; ch0 ch1 ch1; ch1 ch1 ch1};

OPTTH = [];
MAXY = [];
for i1 = 1:8
    color = colors(i1,:);
    plots(end+1) = plot(thresholds,FC(i1,:),'linewidth',lwidth,'color',color);  
    [maxx,maxy] = maxfit(thresholds,FC(i1,:),3);
    line([maxx maxx],[0 0.02],'linewidth',2,'color',color)
    MAXY(end+1) = maxy;
    OPTTH(end+1) = maxx; 
    plot([leg_x0 leg_x0+leg_width],[leg_y0 leg_y0]-(i1-1)*leg_dy,'linewidth',lwidth,'color',color);
    if strcmp(legends{i1,1},ch0)
        text(leg_x0+leg_width+leg_xgap,leg_y0-(i1-1)*leg_dy+leg_yshift+dot_yshift,legends{i1,1},'FontName','Helvetica','fontsize',text_size,'HorizontalAlign','center')
    else
        text(leg_x0+leg_width+leg_xgap,leg_y0-(i1-1)*leg_dy+leg_yshift,legends{i1,1},'FontName','Helvetica','fontsize',text_size,'HorizontalAlign','center')
    end       
    if strcmp(legends{i1,2},ch0)
        text(leg_x0+leg_width+leg_xgap+leg_dx,leg_y0-(i1-1)*leg_dy+leg_yshift+dot_yshift,legends{i1,2},'FontName','Helvetica','fontsize',text_size,'HorizontalAlign','center')
    else
        text(leg_x0+leg_width+leg_xgap+leg_dx,leg_y0-(i1-1)*leg_dy+leg_yshift,legends{i1,2},'FontName','Helvetica','fontsize',text_size,'HorizontalAlign','center')
    end 
    if strcmp(legends{i1,3},ch0)
        text(leg_x0+leg_width+leg_xgap+2*leg_dx,leg_y0-(i1-1)*leg_dy+leg_yshift+dot_yshift,legends{i1,3},'FontName','Helvetica','fontsize',text_size,'HorizontalAlign','center')
    else
        text(leg_x0+leg_width+leg_xgap+2*leg_dx,leg_y0-(i1-1)*leg_dy+leg_yshift,legends{i1,3},'FontName','Helvetica','fontsize',text_size,'HorizontalAlign','center')
    end
end

text(leg_x0+leg_width+leg_xgap+leg_xshift,leg_y0+leg_ygap,'phase','FontName','Helvetica','fontsize',text_size,'rotation',90)
text(leg_x0+leg_width+leg_xgap+leg_dx+leg_xshift,leg_y0+leg_ygap,'spatial period','FontName','Helvetica','fontsize',text_size,'rotation',90)
text(leg_x0+leg_width+leg_xgap+2*leg_dx+leg_xshift,leg_y0+leg_ygap,'contrast','FontName','Helvetica','fontsize',text_size,'rotation',90)

line([thresholds(1) thresholds(end)],[1/K 1/K],'linestyle',':','color','k','linewidth',1)

line([-60 -55],[0 0],'color','k','linewidth',1)

% panel C

axes_position = [mleft+width_A+gapx mbottom width_C height_B];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'FontName','Helvetica','fontsize',text_size,'linewidth',1)

hold on
dx = 0.25;
set(gca,'xlim',[0,8+2*dx],'xticklabels',[])
set(gca,'ylim',[0 1],'ytick',0:.1:1,'yticklabel',[])

MPLEVEL = FC(:,1);
MPLEVEL([2 5 6 8]) = 1/K;
for i1 = 1:8
    patch([i1-1+bargap/2,i1-bargap/2,i1-bargap/2,i1-1+bargap/2]+dx,[0 0 MAXY(i1) MAXY(i1)],colors(i1,:))
    line([i1-1+bargap/2+dx i1-bargap/2+dx],[MPLEVEL(i1) MPLEVEL(i1)],'color','w','linewidth',2)
end
line([1 8]+dx,[1/K 1/K],'linestyle',':','color','k')

xlabel('nuisance combination','FontName','Helvetica','fontsize',text_size)

text(0,1+panel_label_yshift,'C','fontsize',panel_label_size,'FontName','Helvetica','HorizontalAlignment','center')

set(gcf,'PaperPositionMode','auto','papersize',[32 24.5]);
print(gcf,mfilename,'-dpdf','-r0')
saveas(gcf,[mfilename,'.png']);