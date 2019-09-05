% script for creating Fig. 6

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

kappa_values = 1:0.2:2;
noise_values = 0.5:0.5:7;

N_kappa = length(kappa_values);
N_noise = length(noise_values);

OPTRELTH = zeros(N_kappa,N_noise);

kappa_strings = {'1','1d2','1d4','1d6','1d8','2'};
start_index = [3,1,1,1,1,3];
for kappa_index = 1:N_kappa
    load(sprintf('exponent/noise_kappa%s/script_kappa%s_data.mat',kappa_strings{kappa_index},kappa_strings{kappa_index}),'MAXX','V0','V1')
    OPTRELTH(kappa_index,start_index(kappa_index):end) = (MAXX-V0)/V1;
end

font_size = 22;
ABC_size = 34;

ABC = {'A','B','C','D','E'};

color1 = [0 .3 .8];
color2 = [1 .5 0];

width = 232;
height = 262;
mtop = 57;
mbottom = 73;
mleft = 110;
mright = 28;
gapx = 37;

figure_width = 5*width+4*gapx+mleft+mright;
figure_height = height+mbottom+mtop;

set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

% panel A
 
axes_position = [mleft mbottom width height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'linewidth',1,'FontName','Helvetica','fontsize',font_size)

hold on
set(gca,'xlim',[0 0.6],'xtick',0:0.2:1,'ylim',[0 0.75],'ytick',0:0.2:1)

xlabel('noise-to-signal ratio','FontName','Helvetica','fontsize',font_size)
ylab = ylabel({'normalised','optimal threshold'},'FontName','Helvetica','fontsize',font_size)

ylab_pos = get(ylab,'position');
ylab_pos(1) = ylab_pos(1)-0.02;
set(ylab,'position',ylab_pos)

load('V1/phase_nuisance_V1_data.mat','MAXX','V1_values','noise')
optth_values = MAXX;

xx = noise./V1_values;
yy = (optth_values-V0)./V1_values;
scatter(xx,yy,100,'x','MarkerEdgeColor',color1,'linewidth',1.5)

xx2 = noise_values/V1;
yy2 = OPTRELTH(1,:);
scatter(xx2(3:end),yy2(3:end),90,'+','MarkerEdgeColor',color1,'linewidth',1.5)

text(-0.23,0.84,ABC{1},'FontName','Helvetica','fontsize',ABC_size)

% panel BCD

shortname = {'N','K','corr','exponent'};
name = {'N','K','correlation','exponent'};
ref_indices = [3 2 1 1];

line_width = 2;

for n = 1:length(name)-1
    load(sprintf('%s/phase_nuisance_noise3mV_%s_data.mat',shortname{n},shortname{n}))
    eval(sprintf('values = %s_values;',shortname{n}));
    N1 = size(ITERATIONS,1);
    
    axes_position = [mleft+n*width+n*gapx mbottom width height];
    subplot_axes = axes('unit','pixel','position',axes_position);
    set(gca,'linewidth',1,'FontName','Helvetica','fontsize',font_size)
    
    x1 = noise/V1;
    Y1 = (MAXX-V0)/V1;
    
    load(sprintf('%s/phase_nuisance_noise6mV_%s_data.mat',shortname{n},shortname{n}),'MAXX','V0','V1','noise')
    x2 = noise/V1;
    Y2 = (MAXX-V0)/V1;
    
    xshift = 0;
    if n == 3
        xshift = 0.05;
    end
    text(0.07-xshift/2,0.70,name{n},'FontName','Helvetica','fontsize',font_size)
    leg_x = 0.05+xshift/2;
    leg_y = 0.61;
    leg_width = 0.07;
    leg_dy = 0.07;
    leg_dx = 0.013;
    for i = 1:N1
        colori = (N1-i)/(N1-1)*color1+(i-1)/(N1-1)*color2;
        line([x1 x2],[Y1(i) Y2(i)],'color',colori,'linewidth',line_width)
        yline = leg_y-(i-1)*leg_dy;
        line([leg_x leg_x+leg_width],[yline yline],'color',colori,'linewidth',line_width)
        if i == ref_indices(n)
            ref = ' *';
        else
            ref = '';
        end
        text(leg_x+leg_width+leg_dx,yline,sprintf('%g%s',values(i),ref),'FontName','Helvetica','fontsize',font_size)
    end

    set(gca,'xlim',[0 0.6],'xtick',0:0.2:1,'ylim',[0 0.75],'ytick',0:0.2:1,'yticklabel',[])
    xlabel('noise-to-signal ratio','FontName','Helvetica','fontsize',font_size)
    
    text(-0.02,0.84,ABC{n+1},'FontName','Helvetica','fontsize',ABC_size)
end

% panel E

axes_position = [mleft+4*width+4*gapx mbottom width height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'linewidth',1,'fontsize',font_size)

hold on
leg_x = leg_x-0.01;
for i = 1:N_kappa
    colori = (N_kappa-i)/(N_kappa-1)*color1+(i-1)/(N_kappa-1)*color2;
    plot(noise_values(3:end)/V1,OPTRELTH(i,3:end),'linewidth',line_width,'color',colori)
    yline = leg_y-(i-1)*leg_dy;
    line([leg_x leg_x+leg_width],[yline yline],'color',colori,'linewidth',line_width)
    if i == 1
        ref = ' *';
    else
        ref = '';
    end
    text(leg_x+leg_width+leg_dx,yline,sprintf('%g%s',kappa_values(i),ref),'FontName','Helvetica','fontsize',font_size)
end
text(0.07-xshift/2,0.70,name{end},'FontName','Helvetica','fontsize',font_size)

set(gca,'xlim',[0 0.6],'xtick',0:0.2:1,'ylim',[0 0.75],'ytick',0:0.2:1,'yticklabel',[])
xlabel('noise-to-signal ratio','FontName','Helvetica','fontsize',font_size)

text(-0.02,0.84,ABC{5},'FontName','Helvetica','fontsize',ABC_size)

set(gcf,'PaperPositionMode','auto','papersize',[51 14]);
print(gcf,mfilename,'-dpdf','-r0')
saveas(gcf,[mfilename,'.png']);