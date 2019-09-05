% script for creating Fig. 7 figure supplement 1
close all
clearvars

width_A = 320;
width_B = 560;
height = 320;
mtop = 35;
mbottom = 60;
mleft = 65; 
mright = 30;
gapx = 60;
figure_width = mleft+width_A+gapx+width_B+mright;
figure_height = height+mbottom+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

font_size = 16;
ABC_size = 24;

cell_index = 2;
trial_index = 3;
spike_index = 4;

xmin = -4;
xmax = 5;
ymin = -40; 
ymax = 5;

cd(fileparts(which(mfilename)));

% panel A

axes_position = [mleft mbottom width_A height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'linewidth',1,'FontName','Helvetica','fontsize',font_size)

plot([])
hold on

load('data1_raw.mat');

% see search_threshold.m for resource of the following data:
search_thresholds = [2 1.5 1.5 1.5 3 4 1.2];

t_search = 1.5;
n_search = round(t_search/dt);

t1_search = 1.0;
t2_search = 2.0;
n1_search = round(t1_search/dt);
n2_search = round(t2_search/dt);

t3 = 2.5;
t4 = 6;
n3 = round(t3/dt);
n4 = round(t4/dt);

t1_show = 3;
t2_show = 6;
n1_show = round(t1_show/dt);
n2_show = round(t2_show/dt);

fitoptions = optimset('Display','off');

eval([sprintf('U_raw = simple%d_trial%d;',cell_index,trial_index)]);
U_gen = U_raw;
DU = diff(U_raw);
N = length(U_raw);
I_peaks = [-n_search];
SP_thresholds = [];
SP_indices = [];

for i = 1:N-1
    if (DU(i) > search_thresholds(cell_index)) && (i-I_peaks(end) > n_search)
        [peak_value,peak_index] = max(U_raw(i:i+n_search));
        i_peak = i+peak_index-1;
        I_peaks(end+1) = i_peak;

        if (i_peak-n1_search < 1)
            fprintf('simple %d trial %d: start with spike\n',cell_index,trial_index);
            U_gen(1:i_peak+n2_search) = U_gen(i_peak+n2_search);
            SP_thresholds(end+1) = U_gen(i_peak+n2_search);
        elseif (i_peak+n2_search > N-1)
            fprintf('simple %d trial %d: end with spike\n',cell_index,trial_index);
            U_gen(i_peak-n1_search:end) = U_gen(i_peak-n1_search);
            SP_thresholds(end+1) = U_gen(i_peak-n1_search);
        else
            [max_value,max_index] = max(DU(i_peak-n1_search:i_peak));
            if max_index == 1
                fprintf('t1 search warning\n');
            end
            t1 = t1_search-max_index*dt;
            n1 = round(t1/dt);
            t1 = n1*dt;

            [min_value,min_index] = min(DU(i_peak:i_peak+n2_search));
            if min_index == n2_search+1
                fprintf('t2 search warning\n');
            end
            t2 = min_index*dt;
            n2 = min_index;

            factor = 5;
            if n1 < 3
                factor = 7;
            end

            T_fit1 = linspace(-factor*t1,-t1,(factor-1)*n1+1)';
            U_fit1 = U_gen(i_peak-factor*n1:i_peak-n1);
            F1 = @(p,x) (p(2)-(p(1)-x)*p(4)).*(1-sign(x-p(1)))/2 + (p(2)+(x-p(1))*p(3)).*(1+sign(x-p(1)))/2;
            u1 = mean(U_gen(i_peak-factor*n1:i_peak-2*n1));
            hint1 = [-2*t1,u1,(U_gen(i_peak-n1)-u1)/t1,0];
            params1 = lsqcurvefit(F1,hint1,T_fit1,U_fit1,[],[],fitoptions);

            t0 = -params1(1); % ... spike index
            n0 = round(t0/dt);
            spike_threshold = params1(2);
            U_gen(max(i_peak-n0,1):min(i_peak+n2-1,N)) = spike_threshold;

            SP_thresholds(end+1) = spike_threshold;
            SP_indices(end+1) = max(i_peak-n0,1);
            fprintf('.');

            if (i_peak+n4 <= N)
                T_fit2 = linspace(t2,t3,n3-n2+1)';
                U_fit2 = U_raw(i_peak+n2:i_peak+n3);
                u2 = U_raw(i_peak+n2);
                F2 = @(p,x) spike_threshold+(u2-spike_threshold)*exp(-p(1)*(x-t2));
                hint2 = [2.5];
                params2 = lsqcurvefit(F2,hint2,T_fit2,U_fit2,[],[],fitoptions);
                T_remove = linspace(t2,t4,n4-n2+1)';
                U_remove = F2(params2,T_remove)-spike_threshold;
                U_gen(i_peak+n2:i_peak+n4) = U_gen(i_peak+n2:i_peak+n4)-U_remove;
            else
                U_gen(i_peak+n2:end) = spike_threshold;
            end

            if spike_index == length(SP_thresholds)
                n1show = min(2*n1_show,i_peak-1);
                n2show = min(n2_show,N-i_peak);
                T_show = dt*linspace(-n1show,n2show,n1show+n2show+1)';
                U_show = U_raw(i_peak-n1show:i_peak+n2show);
                Y1 = F1(params1,T_fit1);
                Y2 = F2(params2,T_fit2);
                plot(T_show,U_gen(i_peak-n1show:i_peak+n2show),'linewidth',1,'color',[.5 .5 .5])
                plot(T_fit1,Y1,'LineWidth',4,'color',[1 .5 .5])
                plot(T_fit2,Y2,'LineWidth',4,'color',[.5 .5 1])
                scatter(T_fit1(end),U_fit1(end),56,'r')
                scatter(T_fit2(1),U_fit2(1),56,'b')
                scatter(T_show,U_raw(i_peak-n1show:i_peak+n2show),8,'k','filled')
                xlabel('time [ms]','FontName','Helvetica','fontsize',font_size)
                ylabel({'potential [mV]'},'FontName','Helvetica','fontsize',font_size)
                line([-t0 -t0],[-50 spike_threshold],'color','r','linewidth',1);
                line([xmin -t0],[spike_threshold spike_threshold],'color','r','linewidth',1);
                break;    
            end
        end
    end
end

set(gca,'xlim',[xmin xmax],'xtick',xmin:1:xmax,'xticklabel',{'-4','','-2','','0','','2','','4',''})
set(gca,'ylim',[ymin ymax],'ytick',-60:10:20)
                
text(-5.15,ymax+0.05*(ymax-ymin),'A','FontName','Helvetica','fontsize',ABC_size,'HorizontalAlignment','center')

% panel B

axes_position = [mleft+width_A+gapx mbottom width_B height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'linewidth',1,'FontName','Helvetica','fontsize',font_size)

plot([])
hold on

cell_index = 2;
trial_index = 2;

eval([sprintf('U_raw = simple%d_trial%d;',cell_index,trial_index)]);

load('data1_gen.mat');

eval([sprintf('U_gen = simple%d_trial%d_genpot;',cell_index,trial_index)]);

tmin = 23.387;
tmax = 23.454;

n1 = round(1000*tmin/dt);
n2 = round(1000*tmax/dt);

U_raw = U_raw(n1:n2);
U_gen = U_gen(n1:n2);

time = linspace(tmin,tmax,length(U_raw));

plot(time,U_raw,'r','linewidth',2)
plot(time,U_gen,'b','linewidth',2)

set(gca,'xlim',[tmin tmax],'xtick',tmin:0.01:tmax,'xticklabel',{})
set(gca,'ylim',[ymin ymax],'ytick',-60:10:20)

xlabel('time','FontName','Helvetica','fontsize',font_size)

plot([23.437 23.447],[ymin+5 ymin+5],'k','linewidth',2)
text(23.442,ymin+7,'10ms','FontName','Helvetica','fontsize',font_size,'HorizontalAlignment','center')
                
text(tmin-0.003,ymax+0.05*(ymax-ymin),'B','FontName','Helvetica','fontsize',ABC_size,'HorizontalAlignment','center')

set(gcf,'PaperPositionMode','auto','papersize',[24 12])
print(gcf,mfilename,'-dpdf','-r0')
saveas(gcf,[mfilename,'.png']);