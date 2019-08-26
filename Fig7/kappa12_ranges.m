clear all

cd(fileparts(which(mfilename)));

noise_values = 1.5:0.5:7;

N_noise = length(noise_values);

OPTFC = NaN*zeros(2,N_noise);
OPTTH = NaN*zeros(2,N_noise);
OPTRELTH = NaN*zeros(2,N_noise);
THMIN90 = NaN*zeros(2,N_noise);
THMAX90 = NaN*zeros(2,N_noise);
THMIN95 = NaN*zeros(2,N_noise);
THMAX95 = NaN*zeros(2,N_noise);

for kappa = 1:2

    close all

    figure('units','normalized','outerposition',[0 0 1 1],'color','w')
    plot([])
    hold on

    xlabel('FRNL threshold [mV]');
    ylabel('fraction correct');
    title(sprintf('kappa = %d',kappa))

    load(sprintf('script_kappa%d_data.mat',kappa),'thresholds','FC','V0','V1','K')
    
    if kappa == 1
        FC(1,29) = (FC(1,28)+FC(1,30))/2;
    end
    
    if kappa == 2
        FC(1,32) = (FC(1,31)+FC(1,33))/2;
        FC(3,42) = (FC(3,41)+FC(3,43))/2;
    end

    curve_plots = 1:N_noise;
    legend_strings = {};
    for i = 1:N_noise
        color = [1 1-i/N_noise i/N_noise];
        [maxx,maxy,xfit,yfit] = maxfit(thresholds,FC(i,:),4,100);
        line([maxx maxx],[0 0.033],'color',color)
        plot(xfit,yfit,'color',[.8 .8 .8],'linewidth',6);
        curve_plots(i) = plot(thresholds,FC(i,:),'color',color,'linewidth',3);
        legend_strings{end+1} = sprintf('noise = %g mV',noise_values(i));
        %level90 = 0.9*maxy;
        %level95 = 0.95*maxy;
        level90 = 1/K+0.9*(maxy-1/K);
        level95 = 1/K+0.95*(maxy-1/K);
        
        upper90 = find(FC(i,:) > level90);
        upper95 = find(FC(i,:) > level95);
        
        xx = linspace(thresholds(upper90(1)-1),thresholds(upper90(1)),100);
        yy = linspace(FC(i,upper90(1)-1),FC(i,upper90(1)),100);
        xx_upper = xx(yy > level90);
        THMIN90(kappa,i) = xx_upper(1);

        xx = linspace(thresholds(upper90(end)),thresholds(upper90(end)+1),100);
        yy = linspace(FC(i,upper90(end)),FC(i,upper90(end)+1),100);
        xx_upper = xx(yy > level90);
        THMAX90(kappa,i) = xx_upper(end);
        
        xx = linspace(thresholds(upper95(1)-1),thresholds(upper95(1)),100);
        yy = linspace(FC(i,upper95(1)-1),FC(i,upper95(1)),100);
        xx_upper = xx(yy > level95);
        THMIN95(kappa,i) = xx_upper(1);

        xx = linspace(thresholds(upper95(end)),thresholds(upper95(end)+1),100);
        yy = linspace(FC(i,upper95(end)),FC(i,upper95(end)+1),100);
        xx_upper = xx(yy > level95);
        THMAX95(kappa,i) = xx_upper(end);
        
        line([THMIN90(kappa,i) THMAX90(kappa,i)],[level90 level90],'color',color)
        line([THMIN95(kappa,i) THMAX95(kappa,i)],[level95 level95],'color',color,'linestyle',':')
        
        scatter(maxx,maxy,40,'kd','filled')
        scatter(THMIN90(kappa,i),level90,40,'k','filled')
        scatter(THMAX90(kappa,i),level90,40,'k','filled')
        
        OPTFC(kappa,i) = maxy;
        OPTTH(kappa,i) = maxx;
        OPTRELTH(kappa,i) = (maxx-V0)/V1;
            
    end
    legend(curve_plots,legend_strings)
    legend boxoff
 
    saveas(gcf,sprintf('kappa%d_ranges_done.png',kappa));
end

OPTFC(2,end-1) = (OPTFC(2,end-2)+OPTFC(2,end))/2;
OPTTH(2,end-1) = (OPTTH(2,end-2)+OPTTH(2,end))/2;
OPTRELTH(2,end-1) = (OPTTH(2,end-1)-V0)/V1;

save('kappa12_ranges_data.mat','noise_values','V0','V1','OPTFC','OPTTH','OPTRELTH','THMIN90','THMAX90','THMIN95','THMAX95');