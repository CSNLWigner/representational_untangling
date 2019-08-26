% script for creating Fig. 4 figure supplement 2

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

load('script_4d_orientation_data.mat','thresholds','FC')
thresholds0 = thresholds;
FC0_orientation = FC;
load('script_4d_orientation_max_data.mat','thresholds','FC')
thresholds2 = thresholds;
FC2_orientation = mean(FC);

OPTTH_orientation = zeros(1,25);
for i = 1:25
    [maxx,~] = maxfit(thresholds,FC(i,:),3);
    OPTTH_orientation(i) = maxx;
end

thresholds1 = thresholds0(thresholds0 < thresholds2(1));
thresholds3 = thresholds0(thresholds0 > thresholds2(end));
FC1_orientation = FC0_orientation(thresholds0 < thresholds2(1));
FC3_orientation = FC0_orientation(thresholds0 > thresholds2(end));

load('script_4d_wavelength_data.mat','thresholds','FC')
FC0_wavelength = FC;
load('script_4d_wavelength_max_data.mat','thresholds','FC')
FC2_wavelength = mean(FC);

OPTTH_wavelength = zeros(1,25);
for i = 1:25
    [maxx,~] = maxfit(thresholds,FC(i,:),3);
    OPTTH_wavelength(i) = maxx;
end

FC1_wavelength = FC0_wavelength(thresholds0 < thresholds2(1));
FC3_wavelength = FC0_wavelength(thresholds0 > thresholds2(end));

load('script_4d_contrast_data.mat','thresholds','FC')
FC0_contrast = FC;
load('script_4d_contrast_max_data.mat','thresholds','FC')
FC2_contrast = mean(FC);

OPTTH_contrast = zeros(1,25);
for i = 1:25
    [maxx,~] = maxfit(thresholds,FC(i,:),3);
    OPTTH_contrast(i) = maxx;
end

FC1_contrast = FC0_contrast(thresholds0 < thresholds2(1));
FC3_contrast = FC0_contrast(thresholds0 > thresholds2(end));

FC1_contrast(end-2) = (FC1_contrast(end-3)+FC1_contrast(end-1))/2;
FC1_contrast(end-4) = (FC1_contrast(end-5)+FC1_contrast(end-3))/2;

text_size = 16;
ABC_size = 24;
ABC_yshift = 0.08;
lwidth = 3;
ax_width = 1;
xmin = -72.5;
xmax = -41;
ymax = .3;
bargap = 0.22;

orientation_color = [0 .3 .8];
wavelength_color = [0 .7 0]; 
contrast_color = [.9 0 .1];

width_A = 400;
width_B = 100;
height = 320;
mtop = 45;
mbottom = 60;
mright = 27;
mleft = 65;
gapx = 62;

figure_width = mleft+width_A+gapx+width_B+mright;
figure_height = mbottom+height+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

% panel A

axes_position = [mleft mbottom width_A height];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',text_size,'linewidth',ax_width)

hold on
plot(thresholds1,FC1_orientation,'linewidth',lwidth,'color',orientation_color);
plot([thresholds1(end) thresholds2(1)],[FC1_orientation(end) FC2_orientation(1)],'linewidth',lwidth,'color',orientation_color);
plot(thresholds2,FC2_orientation,'linewidth',lwidth,'color',orientation_color);
plot([thresholds2(end) thresholds3(1)],[FC2_orientation(end) FC3_orientation(1)],'linewidth',lwidth,'color',orientation_color);
plot(thresholds3,FC3_orientation,'linewidth',lwidth,'color',orientation_color);  

plot(thresholds1,FC1_wavelength,'linewidth',lwidth,'color',wavelength_color);
plot([thresholds1(end) thresholds2(1)],[FC1_wavelength(end) FC2_wavelength(1)],'linewidth',lwidth,'color',wavelength_color);
plot(thresholds2,FC2_wavelength,'linewidth',lwidth,'color',wavelength_color);
plot([thresholds2(end) thresholds3(1)],[FC2_wavelength(end) FC3_wavelength(1)],'linewidth',lwidth,'color',wavelength_color);
plot(thresholds3,FC3_wavelength,'linewidth',lwidth,'color',wavelength_color);
 
plot(thresholds1,FC1_contrast,'linewidth',lwidth,'color',contrast_color);
plot([thresholds1(end) thresholds2(1)],[FC1_contrast(end) FC2_contrast(1)],'linewidth',lwidth,'color',contrast_color);
plot(thresholds2,FC2_contrast,'linewidth',lwidth,'color',contrast_color);
plot([thresholds2(end) thresholds3(1)],[FC2_contrast(end) FC3_contrast(1)],'linewidth',lwidth,'color',contrast_color);
plot(thresholds3,FC3_contrast,'linewidth',lwidth,'color',contrast_color);

optth_orientation = mean(OPTTH_orientation);
thstd_orientation = std(OPTTH_orientation);
optth_wavelength = mean(OPTTH_wavelength);
thstd_wavelength = std(OPTTH_wavelength);
optth_contrast = mean(OPTTH_contrast);
thstd_contrast = std(OPTTH_contrast);

x1 = optth_orientation-2*thstd_orientation;
x2 = optth_orientation+2*thstd_orientation;
patch([x1 x2 x2 x1],[0 0 .03 .03],orientation_color,'EdgeColor','none')

x1 = optth_wavelength-2*thstd_wavelength;
x2 = optth_wavelength+2*thstd_wavelength;
patch([x1 x2 x2 x1],[0 0 .03 .03],wavelength_color,'EdgeColor','none')

x1 = optth_contrast-2*thstd_contrast;
x2 = optth_contrast+2*thstd_contrast;
patch([x1 x2 x2 x1],[0 0 .03 .03],contrast_color,'EdgeColor','none')
    
K = 10;
line([thresholds1(1) thresholds3(end)],[1/K 1/K],'linestyle',':','linewidth',ax_width,'color','k')

set(gca,'xlim',[xmin,xmax],'xtick',-70:5:-45,'xticklabel',{'-70','-65','-60','-55','-50','-45'})
set(gca,'ylim',[0 ymax],'ytick',0:.1:ymax)

xlabel('FRNL threshold [mV]','FontName','Helvetica','fontsize',text_size)
ylabel('fraction correct','FontName','Helvetica','fontsize',text_size)

x1 = -47;
x2 = -50.5;
y1 = .25;
dy = .025;
dx = 1;
line([x1 x2],[y1 y1]+2*dy,'linewidth',lwidth,'color',orientation_color)
line([x1 x2],[y1 y1]+dy,'linewidth',lwidth,'color',wavelength_color)
line([x1 x2],[y1 y1],'linewidth',lwidth,'color',contrast_color)
text(x1+dx,y1+2*dy,'orientation','FontName','Helvetica','fontsize',text_size)
text(x1+dx,y1+dy,'spatial period','FontName','Helvetica','fontsize',text_size)
text(x1+dx,y1,'contrast','FontName','Helvetica','fontsize',text_size)

text(-75.5,ymax+0.025,'A','fontsize',ABC_size,'FontName','Helvetica','HorizontalAlignment','center','fontweight','bold')

% panel B
 
axes_position = [mleft+width_A+gapx mbottom width_B height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'FontName','Helvetica','fontsize',text_size,'linewidth',1)

hold on
dx = 0.25;
set(gca,'xlim',[0,3+2*dx],'xticklabels',[])
set(gca,'ylim',[0 ymax],'ytick',0:.1:ymax,'yticklabel',[])

[~,maxy_orientation] = maxfit(thresholds2,FC2_orientation,3);
[~,maxy_wavelength] = maxfit(thresholds2,FC2_wavelength,3);
[~,maxy_contrast] = maxfit(thresholds2,FC2_contrast,3);
MAXY = [maxy_orientation,maxy_wavelength,maxy_contrast];
MPLEVEL = [1/K,1/K,1/K];

colors = [orientation_color;wavelength_color;contrast_color];
for i = 1:3
    patch([i-1+bargap/2,i-bargap/2,i-bargap/2,i-1+bargap/2]+dx,[0 0 MAXY(i) MAXY(i)],colors(i,:))
    line([i-1+bargap/2+dx i-bargap/2+dx],[MPLEVEL(i) MPLEVEL(i)],'color','w','linewidth',2)
end
line([0 3]+dx,[1/K 1/K],'linestyle',':','color','k')

xlabel('decoding variable','FontName','Helvetica','fontsize',text_size)

text(0,ymax+0.025,'B','fontsize',ABC_size,'FontName','Helvetica','HorizontalAlignment','center','fontweight','bold')

set(gcf,'PaperPositionMode','auto','papersize',[32 24]);
print(gcf,'Fig4_S2','-dpdf','-r0')
saveas(gcf,sprintf('%s.png',mfilename));