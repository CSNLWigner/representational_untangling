% script for creating Fig. 7 figure supplement 2
% used input data:
%  ../Fig6/exponent/noise_kappa1/script_kappa1_data.mat
%  ../Fig6/exponent/noise_kappa1/script_kappa2_data.mat

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

load('kappa12_ranges_data.mat')

width = 360;
height = 360;
mtop = 47;
mbottom = 60;
mleft = 60; 
mright = 105;
gapx = 35;
figure_width = 2*width+gapx+mleft+mright;
figure_height = height+mbottom+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

font_size = 16;
ABC_size = 26;

for kappa = 1:2

    axes_position = [mleft+(kappa-1)*(width+gapx) mbottom width height];
    subplot_axes = axes('unit','pixel','position',axes_position);
    set(gca,'linewidth',1,'FontName','Helvetica','fontsize',font_size)

    plot([])
    hold on
    load(sprintf('../Fig6/exponent/noise_kappa%d/script_kappa%d_data.mat',kappa,kappa))

    N1 = size(ITERATIONS,1);

    axis([thresholds(1) thresholds(end) 0 1])

    set(gca,'yminortick','on','ytick',0:0.2:1)
    ax = gca;
    ax.YRuler.MinorTickValues = 0:0.1:1;

    xlabel('FRNL threshold [mV]','FontName','Helvetica','fontsize',font_size)

    switch kappa
        case 1
            text(-67.5,1.07,'A','FontName','Helvetica','fontsize',ABC_size)
            text(-60,0.99,'$\kappa = 1$','Interp','Latex','fontsize',ABC_size-6)
            FC(1,29) = (FC(1,28)+FC(1,30))/2-0.012;
        case 2
            text(-70.5,1.07,'B','FontName','Helvetica','fontsize',ABC_size)
            text(-62,0.99,'$\kappa = 2$','Interp','Latex','fontsize',ABC_size-6)
            FC(1,32) = (FC(1,31)+FC(1,33))/2+0.005;
            FC(3,42) = (FC(3,41)+FC(3,43))/2-0.01;
            OPTFC(2,11) = OPTFC(2,11)-0.003;
            THMAX90(2,1) = THMAX90(2,1)+0.065;
    end

    curve_plots = 1:N1;
    legend_strings = {};

    for n1 = 1:N1
        color = [0 1-n1/N1 n1/N1];
        curve_plots(n1) = plot(thresholds,FC(n1,:),'color',color,'linewidth',2);

        line([OPTTH(kappa,n1) OPTTH(kappa,n1)],[0 0.033],'color',color,'linewidth',1.5)
        scatter(OPTTH(kappa,n1),OPTFC(kappa,n1),25,'rd','filled')
        level90 = 1/K+0.9*(OPTFC(kappa,n1)-1/K);
        % line([THMIN90(kappa,n1) THMAX90(kappa,n1)],[level90 level90],'color',color)
        scatter(THMIN90(kappa,n1),level90,15,'k','filled')
        scatter(THMAX90(kappa,n1),level90,15,'k','filled')
        
        legend_strings{end+1} = sprintf('%.2g',noise_values(n1));
    end

    line([thresholds(1) thresholds(end)],[1/K 1/K],'LineStyle',':','Color','k','linewidth',1)

    if kappa == 1
        ylabel('fraction correct','FontName','Helvetica','fontsize',font_size)
    end

    if kappa == 2
        set(gca,'yticklabel',[])
        leg = legend(curve_plots,legend_strings,'FontName','Helvetica','fontsize',font_size);
        legend boxoff
        leg_pos = get(leg,'position');
        leg_pos(1) = leg_pos(1)+0.1;
        set(leg,'position',leg_pos)
        title(leg,'noise [mV]','FontName','Helvetica','fontsize',font_size,'FontWeight','Normal')
        leg.Title.NodeChildren.Position = [0.45 0.98 0];
    end
    
    set(gca,'layer','top')
end

set(gcf,'PaperPositionMode','auto','papersize',[33 17]);
print(gcf,mfilename,'-dpdf','-r0')
saveas(gcf,[mfilename,'.png']);