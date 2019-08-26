% script for creating Fig. 8
% input data: script_fig8_data.mat

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

load('script_fig8_data.mat')

dotsize = 30;
font_size = 22;

width = 800;
height = 600;
mtop = 10;
mbottom = 55;
mright = 5;
mleft = 65;

figure_width = mleft+width+mright;
figure_height = mbottom+height+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

set(gca,'FontName','Helvetica','fontsize',font_size)

maxx = 7;

plot([0 maxx],[PERFORMANCES(1) PERFORMANCES(1)],'linewidth',3,'color',[.75 .75 .75])

hold on

scatter(R_values(2:end),PERFORMANCES(2:end),dotsize,'k','filled')

xfit = R_values(R_values <= maxx);
yfit = PERFORMANCES(R_values <= maxx);

ymin = min(yfit);
ymax = max(yfit);

parabola = @(p,x) p(1)-p(2)*(x.^2);
hint = [ymax,(ymax-ymin)/(maxx^2)];
opt = optimoptions('lsqcurvefit','display','off');
params = lsqcurvefit(parabola,hint,xfit, yfit,[],[],opt);

xx = 0:0.001:maxx;
yy = parabola(params,xx);
plot(xx,yy,'k','linewidth',2)

ymin = 0.8365;
ymax = 0.83925;

set(gca,'ylim',[ymin ymax],'ytick',0.837:0.001:0.839,'fontsize',font_size)
set(gca,'xlim',[0 maxx],'xtick',0:1:7,'fontsize',font_size)

xlabel('threshold perturbation [mV]','FontName','Helvetica','fontsize',font_size+1)
ylabel('fraction correct','FontName','Helvetica','fontsize',font_size+1)

set(gca,'linewidth',2)

hold on

set(gcf,'PaperPositionMode','auto','papersize',[31 25]);
print(gcf,'Fig8','-dpdf','-r0')
saveas(gcf,sprintf('%s.png',mfilename));