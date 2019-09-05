% script for creating Fig. 6 figure supplement 1

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

width = 250;
height = 250;
mtop = 45;
mbottom = 55;
mleft = 60; 
mright = 15;
gapx = 50;
figure_width = 4*width+3*gapx+mleft+mright;
figure_height = height+mbottom+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

ymin = 0.05;
font_size = 16;
ABC_size = 24;

phase_color = [0 .3 .8];
kappa2_color = [1 .5 0];

name = {'N','K','corr','exponent'};
short = {' N ',' K ','    ','    '};
long = {'population size','number of categories','correlation coefficient','FRNL exponent'};
SI_unit = {'' '' '' ''};
colors1 = [102 0 204; 51 102 0; 204 0 102; 0 102 102]/255;
lighting = [1.5 1.5 1.5 1.5];
ymin = [0 0 0 0];
ymax = [1 1 1 1];
ref_ind = [3 2 1 1];
ABCD = {'A','B','C','D'};

cd(fileparts(which(mfilename)));

for i = 1:4

axes_position = [mleft+(i-1)*(width+gapx) mbottom width height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'linewidth',1,'FontName','Helvetica','fontsize',font_size)

plot([])
hold on
if i < 4
    load(sprintf('../Fig6/%s/phase_nuisance_noise3mV_%s_data.mat',name{i},name{i}))
else
    load('phase_nuisance_exponent_data.mat')
end

N1 = size(ITERATIONS,1);

if i == 4
    TEMP = FC_LIN;
    [max_val,max_ind] = max(TEMP(end,:));
    load('phase_nuisance_kappa2_data.mat','FC_LIN')
    TEMP(end,1:max_ind) = FC_LIN(1:max_ind);
    FC_LIN = TEMP;
end

axis([thresholds(1) thresholds(end) ymin(i) ymax(i)])

set(gca,'ytick',0:0.1:1)
set(gca,'yminortick','on','ytick',0:0.1:1)
ax = gca;
ax.YRuler.MinorTickValues = 0:0.05:1;

xlabel('FRNL threshold [mV]','FontName','Helvetica','fontsize',font_size)
if mod(i,4) == 1 ylabel('fraction correct','FontName','Helvetica','fontsize',font_size); end

if i == 1
    text(-46.3,1.013,'$N$','Interp','Latex','fontsize',17)
end
if i == 2
    text(-46.4,1.013,'$K$','Interp','Latex','fontsize',17)
end
if i == 3
    text(-46.1,1.02,'$\varrho$','Interp','Latex','fontsize',18)
end
if i == 4
    text(-45.6,1.02,'$\kappa$','Interp','Latex','fontsize',18)
end

curve_plots = 1:N1;
legend_strings = {};
R = linspace(colors1(i,1),1,lighting(i)*N1);
G = linspace(colors1(i,2),1,lighting(i)*N1);
B = linspace(colors1(i,3),1,lighting(i)*N1);
data(i).optth = [];
data(i).value = [];
for n1 = 1:N1
    if (i == 1) || (i == 2)
        color = [R(N1-n1+1),G(N1-n1+1),B(N1-n1+1)];
    else
        color = [R(n1),G(n1),B(n1)];
    end
    if i == 4
        color = ((N1-n1+1)*phase_color+(n1-1)*kappa2_color)/N1;
    end
    
    curve_plots(n1) = plot(thresholds,FC_LIN(n1,:),'color',color,'linewidth',2);
    if n1 == ref_ind(i) ref = ' *'; else ref = ''; end
    eval(sprintf('value = %s_values(n1);',name{i}));
    legend_strings{end+1} = sprintf('%g%s',value,ref);
    [~,maxi] = max(FC_LIN(n1,:));
    [maxx,maxy] = maxfit(thresholds(maxi-3:maxi+3),FC_LIN(n1,maxi-3:maxi+3));
    line([maxx maxx],[ymin(i) ymin(i)+0.035*(ymax(i)-ymin(i))],'linewidth',2,'color',color)
    data(i).optth(end+1) = maxx;
    data(i).value(end+1) = value;
    if i == 2 line([thresholds(1) thresholds(end)],[1/value 1/value],'LineStyle',':','Color','k','linewidth',1); end
end
if i ~= 2 line([thresholds(1) thresholds(end)],[1/K 1/K],'LineStyle',':','Color','k','linewidth',1); end
leg = legend(curve_plots,legend_strings,'FontName','Helvetica','fontsize',font_size);
leg_pos = get(leg,'position');
leg_pos(1) = leg_pos(1)+0.01;
set(leg,'position',leg_pos)
legend boxoff
text(-42.5,ymax(i)+0.025*(ymax(i)-ymin(i)),[long{i},' (',  '    '   ,') ',SI_unit{i}],'FontName','Helvetica','fontsize',font_size,'horizontalalignment','right');
if i == 1
    text(-77.4,ymin(i)+1.12*(ymax(i)-ymin(i)),ABCD{i},'FontName','Helvetica','fontsize',ABC_size)
else
    text(-75,ymin(i)+1.12*(ymax(i)-ymin(i)),ABCD{i},'FontName','Helvetica','fontsize',ABC_size)
end
end

set(gcf,'PaperPositionMode','auto','papersize',[44 14]);
print(gcf,mfilename,'-dpdf','-r0')
% saveas(gcf,[mfilename,'.png']);