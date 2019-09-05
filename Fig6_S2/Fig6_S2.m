% script for creating Fig. 6 figure supplement 2

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

load('script_2d_ILC_data.mat','thresholds','FC_LIN');
ILC_2D_LIN = FC_LIN;

load('script_2d_Bayes_ILC_data.mat','FC_OPT');
ILC_2D_OPT = FC_OPT;

font_size = 16;
curve_width = 3;
ax_width = 1;

xmin = -72.5;
xmax = -41;

phase_color = [0 .3 .8];
white = [1 1 1];
light_blue = 0.33*phase_color + 0.66*white;
hist_color = [0.66 0.8 0.66];

width = 400;
height = 320;
mtop = 25;
mbottom = 60;
mright = 30;
mleft = 65;

figure_width = mleft+width+mright;
figure_height = mbottom+height+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

axes_position = [mleft mbottom width height];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',ax_width)

plot([])
hold on
 
load('uhist_data.mat','counts','centers');
area(centers,0.45*counts/max(counts),'faceColor',hist_color,'edgecolor','w')    
set(gca,'layer','top')

set(gca,'xlim',[xmin xmax],'xtick',-70:10:-45)
set(gca,'ylim',[0 1],'ytick',0:.1:1,'yticklabel',{'0','','0.2','','0.4','','0.6','','0.8','','1'})

xlabel('FRNL threshold [mV]','FontName','Helvetica','fontsize',font_size)
ylabel('fraction correct','FontName','Helvetica','fontsize',font_size);

p1 = plot(thresholds,ILC_2D_OPT,'linewidth',curve_width,'color',light_blue);
p2 = plot(thresholds,ILC_2D_LIN,'linewidth',curve_width,'color',phase_color);

[maxx,maxy] = maxfit(thresholds,ILC_2D_LIN);
scatter(maxx,maxy,100,phase_color,'d','filled')

leg = legend([p1 p2],{' optimal decoder',' linear decoder'},'location','northwest');
a = get(leg,'position');
a(1) = a(1)+0.05;
set(leg,'position',a);
legend boxoff

gwidth = 69;
axes_position = [mleft+width-1.65*gwidth mbottom+height-1.1*gwidth gwidth gwidth];
subplot_axes = axes('unit','pixel','position',axes_position);
hold on
axis off
daspect([1 1 1])

gabor_params = [0,0,3,45,0,10];
NN = 1000;
G = circular_gabor([5],gabor_params,NN);
maxg = max(max(G));
for i = 1:NN
    for j = 1:NN
        if (i-NN/2)^2+(j-NN/2)^2 > (NN/2)^2
            G(i,j) = maxg;
        end
    end 
end  
imagesc(G)
colormap(0.25+0.75*gray(256));

cx0 = NN/2;
cy0 = NN/2;
rr = 1.5*NN/2;
ph0 = pi/4;
dph = pi*12/180;
tt = linspace(ph0-dph,ph0+dph);
xx = cx0+rr*cos(tt);
yy = cy0+rr*sin(tt);
hold on
plot(xx,yy,'color','k','linewidth',1.5)

asize = 90;
tilt1 = 20*pi/180;
tilt2 = 25*pi/180;
ax0 = cx0+rr*cos(ph0-dph);
ay0 = cy0+rr*sin(ph0-dph);
plot([ax0 ax0-asize*cos(pi/2-ph0+dph+tilt1)],[ay0 ay0+asize*sin(pi/2-ph0+dph+tilt1)],'color','k','linewidth',1.5)
plot([ax0 ax0-asize*cos(pi/2-ph0+dph-tilt2)],[ay0 ay0+asize*sin(pi/2-ph0+dph-tilt2)],'color','k','linewidth',1.5)
ax0 = cx0+rr*cos(ph0+dph);
ay0 = cy0+rr*sin(ph0+dph);
plot([ax0 ax0+asize*cos(pi/2-ph0-dph+tilt2)],[ay0 ay0-asize*sin(pi/2-ph0-dph+tilt2)],'color','k','linewidth',1.5)
plot([ax0 ax0+asize*cos(pi/2-ph0-dph-tilt1)],[ay0 ay0-asize*sin(pi/2-ph0-dph-tilt1)],'color','k','linewidth',1.5)

text(1.44*rr,1.56*rr,sprintf('10%c',char(176)),'FontName','Helvetica','fontsize',font_size-2)

set(gcf,'PaperPositionMode','auto','papersize',[18 16]);
print(gcf,mfilename,'-dpdf','-r0')
saveas(gcf,[mfilename,'.png']);