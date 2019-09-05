% script for creating Fig. 3 figure supplement 3
% extrenal functions used: -
% input data:
%   script_2d_methods_data.mat (phase nuisance, different methods)
%   script_2d_Bayes_data.mat (phase nuisance, Bayes)

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

load('script_2d_methods_data.mat','thresholds','FC_LIN','PFC_LIN','COS_LIN');

load('../Fig3/script_2d_Bayes_data.mat','FC_OPT','PFC_OPT','K');
MI = log(PFC_OPT)+log(K);

text_size = 16;
ABC_size = 24;
lwidth = 3;

xmin = -72.5;
xmax = -41;

width_A = 480;
height = 320;
gapx = 155;
mtop = 55;
mbottom = 62;
mright = 75;
mleft = 95;

figure_width = mleft+2*width_A+gapx+mright;
figure_height = mbottom+height+mtop;

fig = figure;
black = [0 0 0];
orange = [1 .55 0];
set(fig,'unit','pixel','position',[0 0 figure_width figure_height])
set(fig,'color','white','defaultAxesColorOrder',[black;orange])

axes_position = [mleft mbottom width_A height];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',text_size,'linewidth',1)

axis([xmin xmax 0 1])
xlabel('FRNL threshold [mV]');

yyaxis left
hold on
plot(thresholds,FC_LIN,'color','b','linewidth',2,'linestyle','-')
plot(thresholds,PFC_LIN,'color','k','linewidth',2,'linestyle','-')
line([thresholds(1) thresholds(end)],[1/K 1/K],'color','r','linestyle',':');

ylabel({'\color{blue}fraction correct','\color{black}probabilistic fraction correct'});

yyaxis right

plot(thresholds,0.1+0.9*COS_LIN,'linewidth',2,'color',orange);

ylabel('\color[rgb]{1 .55 0}1- cosine error');

yticks(.1:.9/5:1)
yticklabels({'0','0.2','0.4','0.6','0.8','1'})

[maxx_fc,maxy_fc] = maxfit(thresholds,FC_LIN,3);
[maxx_pfc,maxy_pfc] = maxfit(thresholds,PFC_LIN,3);
[maxx_cos,maxy_cos] = maxfit(thresholds,COS_LIN,3);
line([maxx_fc maxx_fc],[0 0.033],'color','b','linewidth',2)
line([maxx_pfc maxx_pfc],[0 0.033],'color','k','linewidth',2)
line([maxx_cos maxx_cos],[0 0.033],'color',orange,'linewidth',2)
scatter(maxx_fc,maxy_fc,100,'r','d','filled')
scatter(maxx_pfc,maxy_pfc,100,'r','d','filled')
scatter(maxx_cos,0.1+0.9*maxy_cos,100,'r','d','filled')

title({'linear decoder',''})

text(-76.5,1+0.11,'A','fontsize',ABC_size,'HorizontalAlignment','center')

axes_position = [mleft+width_A+gapx mbottom width_A height];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',text_size,'linewidth',1)

xlabel('FRNL threshold [mV]');

yyaxis left
hold on
plot(thresholds,FC_OPT,'color','b','linewidth',2,'linestyle','-')
plot(thresholds,PFC_OPT,'color','k','linewidth',3,'linestyle','-')
line([thresholds(1) thresholds(end)],[1/K 1/K],'color','r','linestyle',':');

axis([xmin xmax 0 1])

ylabel({'\color{blue}fraction correct','\color{black}probabilistic fraction correct'});

yyaxis right

hold on
plot(thresholds,exp(MI),'linewidth',1.5,'color',orange);

axis([xmin xmax 0 K])

ylabel('\color[rgb]{1 .55 0}mutual information [nat]');

yticks(exp(0:.5:2))
yticklabels({'0','0.5','1','1.5','2'})

minorticks = exp(0:.1:2.3);
for i = 1:length(minorticks)
    line([-41.2 -41],[minorticks(i) minorticks(i)],'color',orange,'linewidth',1)
end
     
title({'optimal decoder',''})

text(-76.5,11.1,'B','fontsize',ABC_size,'HorizontalAlignment','center')

set(gcf,'PaperPositionMode','auto','papersize',[45 16]);
print(gcf,mfilename,'-dpdf','-r0')
saveas(gcf,[mfilename,'.png']);