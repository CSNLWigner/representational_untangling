load('data_raw.mat');

show = false;
save_figure = false;

% see search_threshold.m for resource of the following data:
search_thresholds = [4 2 2 2 3 4 1.2];

cd(fileparts(which(mfilename)));
mdatafile = 'data_gen.mat';
save(mdatafile,'dt');

[t_search t1 t2 t3 t4 t5] = deal(1.5,0.8,0.8,0.45,2.5,6);
n_values = num2cell(round([t_search t1 t2 t3 t4 t5]/dt));
[n_search n1 n2 n3 n4 n5] = deal(n_values{:});

figure('units','normalized','outerposition',[0 0 1 1])

for n = 1:size(cells,2)
    m = sum(strcmp({cells(1:n).type},cells(n).type));
    for trial = 1:cells(n).trials
        load(sprintf('data/data_%s_00%d.mat',cells(n).id,trial));
        eval([sprintf('U_raw = %s%d_trial%d;',cells(n).type,m,trial)]);
        U_gen = U_raw;
        DU = diff(U_raw);
        N = length(U_raw);
        DU_peaks = findpeaks(DU);
        hist(DU_peaks,100)
        
       
        break
        
        I_peaks = [-n_search];
        SP_thresholds = [];
        for i = 1:N-1
            if (DU(i) > search_thresholds(n)) && (i-I_peaks(end) > n_search)
                [peak_value,peak_index] = max(U_raw(i:i+n_search));
                i_peak = i+peak_index-1;
                I_peaks(end+1) = i_peak;
                if (i_peak <= n1)
                    fprintf('%s %d trial %d: start with spike',cells(n).type,m,trial);
                    SP_thresholds(end+1) = NaN;
                elseif (i_peak > N-n2)
                    fprintf('%s %d trial %d: end with spike',cells(n).type,m,trial);
                    U_gen(i_peak-n1:end) = U_gen(i_peak-n1);
                    SP_thresholds(end+1) = U_gen(i_peak-n1);
                else               
                    T_fit1 = linspace(-t1,t2,n1+n2+1)';
                    U_fit1 = U_raw(i_peak-n1:i_peak+n2);
                    F1 = @(p,x) p(6)+p(5)*(1+sign(p(2)*x-p(1))).*exp(-(log(p(2)*x-p(1))-p(3)).^2/p(4)/p(4)/2)/sqrt(2*pi)/p(4)./(p(2)*x-p(1));
                    u_background = U_raw(i_peak-n1);
                    hint1 = [-0.9206,0.6266,0.1661,0.5614,24.8021,u_background];           
                    params1 = lsqcurvefit(F1,hint1,T_fit1,U_fit1);
                    t0 = -params1(1)/params1(2);
                    n0 = round(t0/dt);
                    if (i_peak-n0 > 0)
                        spike_threshold = U_raw(i_peak-n0);
                        U_gen(i_peak-n0:min(i_peak+n3-1,N)) = spike_threshold;
                    else
                        spike_threshold = U_gen(i_peak-n1);
                        U_gen(1:i_peak+n3-1) = U_gen(1);
                    end 
                    SP_thresholds(end+1) = spike_threshold;
                    if (i_peak+n5 <= N)
                        T_fit2 = linspace(t3,t4,n4-n3+1)';
                        U_fit2 = U_raw(i_peak+n3:i_peak+n4);
                        u3 = U_raw(i_peak+n3);
                        F2 = @(p,x) spike_threshold+(u3-spike_threshold)*exp(-p(1)*(x-t3));
                        hint2 = [2.5];
                        params2 = lsqcurvefit(F2,hint2,T_fit2,U_fit2);
                        T_remove = linspace(t3,t5,n5-n3+1)';
                        U_remove = F2(params2,T_remove)-spike_threshold;
                        U_gen(i_peak+n3:i_peak+n5) = U_gen(i_peak+n3:i_peak+n5)-U_remove;
                    else
                        U_gen(i_peak+n3:end) = spike_threshold;
                    end
                    if show
                        n_show = min(2*n1,i_peak-1);
                        n_remove = min(n5,N-i_peak);
                        T_show = dt*linspace(-n_show,n_remove,n_show+n_remove+1)';
                        U_show = U_raw(i_peak-n_show:i_peak+n_remove);
                        Y1 = F1(params1,T_fit1);
                        Y2 = F2(params2,T_fit2);
                        plot(T_show,U_gen(i_peak-n_show:i_peak+n_remove),'k')
                        hold on
                        rectangle('Position',[-t1 spike_threshold-2 t1+t2 1],'FaceColor','r');     
                        rectangle('Position',[t3 spike_threshold-3 t4-t3 0.6],'FaceColor','b');
                        plot(T_fit1,Y1,'r','LineWidth',5)
                        plot(T_fit2,Y2,'b','LineWidth',3)
                        plot(T_show,U_show,'k.')
                        % title(sprintf('spike %d (peak index %d/%d)',length(I_peaks)-1,i_peak,N))
                        xlabel('time [ms]')
                        ylabel('potential [mV]')
                        line([-t0 -t0],[spike_threshold-3 spike_threshold],'color','k');
                        line([0 0],[spike_threshold-3 spike_threshold],'color','k');
                        %saveas(gcf,'spike.pdf');
                        %print('spike','-dpng','-r300');
                        hold off
                        w = waitforbuttonpress;
                        if w == 0
                            break
                        end
                    end
                end
            end
        end
        varname_data = sprintf('%s%d_trial%d_peaks',cells(n).type,m,trial);
        eval([sprintf('%s = I_peaks(2:end);',varname_data)]);
        save(mdatafile,varname_data,'-append');
        varname_data = sprintf('%s%d_trial%d_thresholds',cells(n).type,m,trial);
        eval([sprintf('%s = SP_thresholds;',varname_data)]);
        save(mdatafile,varname_data,'-append');
        varname_data = sprintf('%s%d_trial%d_genpot',cells(n).type,m,trial);
        eval([sprintf('%s = U_gen;',varname_data)]);
        save(mdatafile,varname_data,'-append');
    end
    
    break
    
end