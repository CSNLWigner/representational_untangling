% script for creating Fig. 7

close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

font_size = 20;
ABC_size = 30;
math_size = 21;

yellow = [255 215 0]/255;
red = [1 .3 .3];
white = [1 1 1];
cell_color = [30 144 255]/255;
grey = [.84 .84 .84];

mpdistColor = [127 127 127]/255;
spikeColor = 0.66*red+0.33*yellow;

legendColor = [1 0 .0];
light_orange = 0.33*mpdistColor+0.66*white;
firingrateColor = [0 0 0];
firingratefitColor = [.5 .5 .5];

width_A = 280;
width_B = 250;
gapx1 = 170;
height_up = 150;
height_dn = 80;
gapy1 = 35;
gapy2plus = 20;
height_B = height_dn+gapy1+height_up+gapy2plus;
gapx2 = 105;
width_C = (width_A+gapx1+width_B-gapx2)/2;
width_Cplus = 25;
height_C = 250;
gapy2 = 120;
gapx3 = 40;
width_legend = 225;
legend_yshift = 20;
mtop = 58;
mbottom = 70;
mright = 0;
mleft = 87;
figure_width = mleft+width_A+gapx1+width_B+gapx3+width_legend+mright;
figure_height = mbottom+height_B+gapy2+height_C+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

% panel A

axes_position = [mleft mbottom+height_C+gapy2 width_A height_B];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

plot([])
hold on

% load raw membrane potential data for illustration
load('MP_data_for_illustration.mat','u_raw','u_gen')

NC = 6; % number of cycles
TC = 0.5; % time duration of a single cycle
ST = NC*TC; % time duration of the stimulus
L = length(u_raw);
dt = ST/(L-1);
time = 0:dt:ST;
k = 1.25;
t1 = 0.325;
n1 = round(t1/dt)+1;
nk = round(k*TC/dt);
n2 = n1+nk;

A = 8.5;
xx = time(n1:n2);
DC = mean(u_gen(n1:n2));
yy = DC+A*sin(2*pi*xx/TC+1.87);

error = u_gen(n1:n2)'-yy;
lower_error = error(yy < DC);
error_std = std(lower_error);

fill([time(n1:n2),fliplr(time(n1:n2))],[yy+error_std,fliplr(yy-error_std)],light_orange,'edgecolor','none')

plot(xx,yy,'color',mpdistColor,'linewidth',3)
plot([time(n1) time(n2)],[DC DC],'linestyle','-','color',legendColor,'linewidth',2)

plot(time(n1:n2),u_raw(n1:n2),'color',spikeColor)
plot(time(n1:n2),u_gen(n1:n2),'k','linewidth',2)

set(gca,'linewidth',1,'FontName','Helvetica','fontsize',font_size)
set(gca,'xcolor','w')
set(gca,'xlim',[t1-0.025 t1+k*TC])
set(gca,'ylim',[-57 -30],'ytick',-100:10:0)

ylab = ylabel('membrane potential [mV]','FontName','Helvetica','fontsize',font_size);
set(ylab, 'position', get(ylab,'position')-[0.01,0,0]);

threshold = -40.2429; % [mV]

plot([time(n1) time(n2)],[threshold threshold],'linestyle','-','color',spikeColor,'linewidth',2)

tunit = 100/1000;
thalf = (t1+t1+k*TC)/2;
shift = -180/1000;
line([thalf-tunit/2 thalf+tunit/2]+shift,[min(u_gen) min(u_gen)],'linewidth',5,'color','k')
text(thalf+shift-0.038,min(u_gen)+2.0,'100','Fontname','Helvetica','fontsize',math_size-3,'horizontalalignment','center')
text(thalf+shift+0.040,min(u_gen)+2.0,'ms','Fontname','Helvetica','fontsize',math_size-3,'horizontalalignment','center')

arrow_x = t1+k*TC+0.025
text(arrow_x+0.035,DC,'2','FontName','Helvetica','fontsize',math_size,'horizontalalignment','center','color',legendColor)
text(arrow_x+0.065,DC,'u','FontName','Helvetica','fontsize',math_size,'horizontalalignment','center','color',legendColor)
text(arrow_x+0.105,DC-0.46,'AC','FontName','Helvetica','fontsize',math_size-6,'horizontalalignment','center','color',legendColor)

arrow_x = thalf+0.09;
arrow_y = DC-A;
text(arrow_x+0.044,arrow_y-1.75*error_std,'$$\pm$$','FontName','Helvetica','fontsize',math_size,'interpreter','LaTex','horizontalalignment','center','color',legendColor)
text(arrow_x+0.088,arrow_y-1.71*error_std,'$$\sigma$$','FontName','Helvetica','fontsize',math_size+5,'interpreter','LaTex','horizontalalignment','center','color',legendColor)

label_x = thalf+0.025;
yshift = 1.55;

text(label_x+0.045,threshold+yshift,'u','FontName','Helvetica','fontsize',math_size,'horizontalalignment','center','color',spikeColor)
text(label_x+0.077,threshold+yshift-0.46,'th','FontName','Helvetica','fontsize',math_size-6,'horizontalalignment','center','color',spikeColor)

text(label_x+0.045,DC+yshift,'u','FontName','Helvetica','fontsize',math_size,'horizontalalignment','center','color',legendColor)
text(label_x+0.087,DC+yshift-0.47,'DC','FontName','Helvetica','fontsize',math_size-6,'horizontalalignment','center','color',legendColor)

text(-0.233,1.125,'A','fontname','Helvetica','fontsize',ABC_size,'units','normalized')

arrow([.359 .71],[.359 .83],.015,.011,2.2,'r')
arrow([.359 .71],[.359 .60],.015,.011,2.2,'r')

arrow([.257 .543],[.257 .578],.015,.011,2.2,'r')
arrow([.257 .661],[.257 .624],.015,.011,2.2,'r')

% panel B dn (histograms)

axes_position = [mleft+width_A+gapx1 mbottom+height_C+gapy2+gapy2plus width_B height_dn];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

tbin_ms = 10;
thin_line = 1.5;
thick_line = 3;
dt = 8.0550e-05;

% load spike counting data for illustration
load('SP_data_for_illustration.mat','dt','umin','umax','du','UBIN_CENTERS','tbin','nbin','U_HIST','SP1_HIST')

plot([])
hold on

normhist_U = U_HIST/sum(U_HIST);
normhist_SP1 = SP1_HIST/sum(SP1_HIST);

xx = [UBIN_CENTERS-du/2,UBIN_CENTERS(end)+du/2]; xx = [xx;xx]; xx = xx(:);

yy = [normhist_U;normhist_U]; yy = yy(:); yy = [0;yy;0];
p1 = area(xx,yy,'faceColor',mpdistColor,'edgecolor',mpdistColor,'linewidth',2);
set(p1,'facealpha',.2) 

yy_sp1 = [normhist_SP1;normhist_SP1]; yy_sp1 = yy_sp1(:); yy_sp1 = [0;yy_sp1;0];
p2 = area(xx(40:62),yy_sp1(40:62),'faceColor',spikeColor,'edgecolor',spikeColor,'linewidth',2);
set(p2,'facealpha',.2)

xmin = -62;
xmax = -33;
xlim([xmin xmax]);
ymax = 1.1*max(max(normhist_U),max(normhist_SP1));
ylim([0 ymax]);
set(gca,'ytick',[])

xlabel('membrane potential [mV]','FontName','Helvetica','fontsize',font_size)
text(-68.15,-0.01,'normalised','FontName','Helvetica','fontsize',font_size,'rotation',90)
text(-65.25,0.045,'counts','FontName','Helvetica','fontsize',font_size,'rotation',90)

threshold = -44.8730;

plot([threshold threshold],[0 ymax],'color',spikeColor,'linewidth',2,'linestyle','-');

% panel gap

axes_position = [mleft+width_A+gapx1 mbottom+height_C+gapy2+height_dn+gapy2plus width_B gapy1];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

plot([threshold threshold],[0 1],'color',spikeColor,'linewidth',2,'linestyle','-');

xmin = -62;
xmax = -33;
xlim([xmin xmax]);

ylim([0 1])

axis off

% panel B up

axes_position = [mleft+width_A+gapx1 mbottom+height_C+gapy2+height_dn+gapy1+gapy2plus width_B height_up];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

plot([])
hold on

FRNL = SP1_HIST./U_HIST/tbin;
FRNL_max = max(FRNL);

FRNL(U_HIST<30) = NaN;

FRNL_ERR = sqrt(SP1_HIST.*(U_HIST-SP1_HIST)./U_HIST)./U_HIST/tbin;
FRNL_ERR(U_HIST<30) = NaN;
FRNL_ERR = 2*FRNL_ERR;

kappa_values = [1 2];
prefactor_values = [];
for k = 1:length(kappa_values)
    unnorm_R = (UBIN_CENTERS-threshold).^kappa_values(k);
    unnorm_R(UBIN_CENTERS < threshold) = 0;
    prefactor_values(end+1) = sum(SP1_HIST)/tbin/sum(U_HIST.*unnorm_R); 
end
 
xmin = -62;
xmax = -33;
xlim([xmin xmax]);
set(gca,'xticklabel',[])

ymax = 100;
ylim([0 ymax]);
set(gca,'ytick',0:25:ymax)

ylabel(' firing rate [Hz]','FontName','Helvetica','fontsize',font_size)

F1 = FRNL-FRNL_ERR;
F2 = FRNL+FRNL_ERR;
F1 = F1(6:28);
F2 = F2(6:28);
UU = UBIN_CENTERS(6:28);
p3 = fill([UU,fliplr(UU)],[F1,fliplr(F2)],[.75 .75 .75],'edgecolor','none');
set(p3,'facealpha',.5)

xxx = linspace(xmin,xmax,100);
for k = 1:length(kappa_values)
    kappa = kappa_values(k);
    prefactor = prefactor_values(k);
    yyy = prefactor*(xxx-threshold).^kappa;
    yyy(xxx < threshold) = 0;
    plot(xxx,yyy,'color',legendColor,'linewidth',thin_line,'linestyle','-')
end

plot(UBIN_CENTERS,FRNL,'color',firingrateColor,'linewidth',thick_line,'linestyle','-')

plot([threshold threshold],[0 ymax],'color',spikeColor,'linewidth',2,'linestyle','-');

text(threshold-0.8,ymax+11,'u','FontName','Helvetica','fontsize',math_size,'horizontalalignment','center','color',spikeColor)
text(threshold+0.85,ymax+8,'th','FontName','Helvetica','fontsize',math_size-6,'horizontalalignment','center','color',spikeColor)

text(-35,ymax+10,'$$\kappa_1$$','FontName','Helvetica','fontsize',math_size,'interpreter','LaTex','horizontalalignment','center','color',legendColor)
text(-38.5,ymax+10,'$$\kappa_2$$','FontName','Helvetica','fontsize',math_size,'interpreter','LaTex','horizontalalignment','center','color',legendColor)

text(-0.345,1.22,'B','fontname','Helvetica','fontsize',ABC_size,'units','normalized')

% panel C

axes_position = [mleft mbottom width_C+width_Cplus height_C];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

plot([])
hold on

% load neural data fitting results
load('c100_segments.mat');

xrange = [0 0.6];
yrange = [0 0.95];

load('results.mat');

DC_c100 = [-46.7742556175103;-42.3531879290987;-47.2701417902820;-36.9291760170919];

DC_values = DC_c100;
AC_values = zeros(1,4);
noise_values = zeros(1,4);
nonlinth_values = zeros(1,4);

for n = 1:4
    cell_AC = [];
    cell_noise = [];
    for i = 1:length(segment)
        if segment(i).cell_index == n  
            cell_AC = [cell_AC segment(i).sinfit_AC];
            cell_noise = [cell_noise segment(i).lowernoise_20ms];
            cell_nonlinth = segment(i).nonlinth;
        end
    end
    AC_values(n) = max(cell_AC);
    noise_values(n) = mean(cell_noise);
    nonlinth_values(n) = cell_nonlinth;
end

hold on
set(gca,'FontName','Helvetica','fontsize',font_size)
set(gca,'xlim',xrange,'ylim',yrange,'fontsize',font_size,'linewidth',1.5)

title({'representational untangling',''})

V0 = -60; V1 = 12;

sigma_values = 1.5:0.5:7;

load('kappa12_ranges_data.mat','THMIN90','THMAX90')

y1 = THMIN90(1,:);
y2 = THMAX90(2,:);
fill([sigma_values,fliplr(sigma_values)]/V1,[y1-V0,fliplr(y2)-V0]/V1,grey,'edgecolor','none')
set(gca,'layer','top') 

for n = 1:4
    x = noise_values(n)/AC_values(n);
    y = (nonlinth_values(n)-DC_values(n))/AC_values(n);
    scatter(x,y,200,cell_color,'filled')
end

xlabel('noise-to-signal ratio','FontName','Helvetica','fontsize',font_size)
ylab = ylabel('normalised threshold','FontName','Helvetica','fontsize',font_size)
set(ylab, 'position', get(ylab,'position')-[0.01,0,0]);

text(-0.205,1.175,'C','fontname','Helvetica','fontsize',ABC_size,'units','normalized')

axes_position = [mleft+0.15*width_C mbottom+0.7*height_C 0.4*width_C 0.3*height_C];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

hold on
fill([sigma_values,fliplr(sigma_values)]/V1,[y1-V0,fliplr(y2)-V0]/V1,grey,'edgecolor','none')
set(gca,'xlim',[0 .6],'xtick',0:.1:.6,'xticklabel',[])
set(gca,'ylim',[-2 3],'ytick',-2:1:3,'yticklabel',[])
set(gca,'layer','top')

text(-0.03,-2,'-2','fontname','Helvetica','fontsize',13,'HorizontalAlignment','right')
text(-0.03,0,'0','fontname','Helvetica','fontsize',13,'HorizontalAlignment','right')
text(-0.03,3,'3','fontname','Helvetica','fontsize',13,'HorizontalAlignment','right')
text(0.1,-2.6,'0.1','fontname','Helvetica','fontsize',13,'HorizontalAlignment','center')
text(0.6,-2.6,'0.6','fontname','Helvetica','fontsize',13,'HorizontalAlignment','center')

all_points = 0;
grey_region = 0;

for i1 = 1:4
DC = DC_values(i1);
for i2 = 1:4
AC = AC_values(i2);
for i3 = 1:4
noise = noise_values(i3);
for i4 = 1:4
nonlinth = nonlinth_values(i4);
x = noise/AC;
y = (nonlinth-DC)/AC;
if ~((i1 == i2) && (i2 == i3) && (i3 == i4))
    if x < 0.6
        scatter(x,y,5,'k','filled')
        all_points = all_points+1;
        [~,closest] = min(abs(x-sigma_values/V1));
        if ((y1(closest)-V0)/V1 <= y) && (y <= (y2(closest)-V0)/V1)
           grey_region = grey_region+1;
        end
    end
end
end
end
end
end

% panel D

load('script_2d_Bayes_MP_data','noise_values','FC_OPT','K')
Bayes_noise_values = noise_values;
Bayes_MPFC = FC_OPT;
Bayes_MPFC(1:3) = 1;
Bayes_MPFC90 = 1/K+.9*(Bayes_MPFC-1/K);
TH90 = NaN*noise_values;
noise_strings = {'1d8','2d4','3d0','3d6','4d2','4d8','5d4','6d0','6d6','7d2','7d8'};
for i = 1:length(noise_values)
    load(sprintf('script_2d_Bayes_noise%smV_data.mat',noise_strings{i}),'thresholds','FC_OPT')
    
    fprintf('noise = %gmV, OPT FC MAX = %g, level90 = %g, scan range: [%g,%g]\n',noise_values(i),Bayes_MPFC(i),Bayes_MPFC90(i),FC_OPT(1),FC_OPT(end))
    
    for j = 1:length(thresholds)
        if FC_OPT(j) < Bayes_MPFC90(i)
            threshold1 = thresholds(j-1);
            threshold2 = thresholds(j);
            fc1 = FC_OPT(j-1);
            fc2 = FC_OPT(j);
            TH90(i) = threshold1+(threshold2-threshold1)*(fc1-Bayes_MPFC90(i))/(fc1-fc2);
            break
        end
    end
end

axes_position = [mleft+width_C+gapx2 mbottom width_C+width_Cplus height_C];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

plot([])
hold on

% load neural data fitting results
load('c100_segments.mat');

xrange = [0 0.6];
yrange = [0 0.95];

load('results.mat');

DC_c100 = [-46.7742556175103;-42.3531879290987;-47.2701417902820;-36.9291760170919];

DC_values = DC_c100;
AC_values = zeros(1,4);
noise_values = zeros(1,4);
nonlinth_values = zeros(1,4);

for n = 1:4
    cell_AC = [];
    cell_noise = [];
    for i = 1:length(segment)
        if segment(i).cell_index == n  
            cell_AC = [cell_AC segment(i).sinfit_AC];
            cell_noise = [cell_noise segment(i).lowernoise_20ms];
            cell_nonlinth = segment(i).nonlinth;
        end
    end
    AC_values(n) = max(cell_AC);
    noise_values(n) = mean(cell_noise);
    nonlinth_values(n) = cell_nonlinth;
end

hold on
set(gca,'FontName','Helvetica','fontsize',font_size)
set(gca,'xlim',xrange,'ylim',yrange,'fontsize',font_size,'linewidth',1.5)

title({'information maximisation',''})

V0 = -60; V1 = 12;

sigma_values = 1.5:0.5:7;

load('kappa12_ranges_data.mat','THMIN90','THMAX90')

y1 = 0*Bayes_noise_values-100;
y2 = TH90;
fill([Bayes_noise_values,fliplr(Bayes_noise_values)]/V1,[y1-V0,fliplr(y2)-V0]/V1,grey,'edgecolor','none')
set(gca,'layer','top') 

for n = 1:4
    x = noise_values(n)/AC_values(n);
    y = (nonlinth_values(n)-DC_values(n))/AC_values(n);
    scatter(x,y,200,cell_color,'filled')
end

xlabel('noise-to-signal ratio','FontName','Helvetica','fontsize',font_size)

text(-0.113,1.175,'D','fontname','Helvetica','fontsize',ABC_size,'units','normalized')

axes_position = [mleft+width_C+gapx2+0.15*width_C mbottom+0.7*height_C 0.4*width_C 0.3*height_C];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)

hold on
fill([Bayes_noise_values,fliplr(Bayes_noise_values)]/V1,[y1-V0,fliplr(y2)-V0]/V1,grey,'edgecolor','none')
set(gca,'xlim',[0 .6],'xtick',0:.1:.6,'xticklabel',[])
set(gca,'ylim',[-2 3],'ytick',-2:1:3,'yticklabel',[])
set(gca,'layer','top')

text(-0.03,-2,'-2','fontname','Helvetica','fontsize',13,'HorizontalAlignment','right')
text(-0.03,0,'0','fontname','Helvetica','fontsize',13,'HorizontalAlignment','right')
text(-0.03,3,'3','fontname','Helvetica','fontsize',13,'HorizontalAlignment','right')
text(0.1,-2.6,'0.1','fontname','Helvetica','fontsize',13,'HorizontalAlignment','center')
text(0.6,-2.6,'0.6','fontname','Helvetica','fontsize',13,'HorizontalAlignment','center')

all_points = 0;
grey_region = 0;

for i1 = 1:4
DC = DC_values(i1);
for i2 = 1:4
AC = AC_values(i2);
for i3 = 1:4
noise = noise_values(i3);
for i4 = 1:4
nonlinth = nonlinth_values(i4);
x = noise/AC;
y = (nonlinth-DC)/AC;
if ~((i1 == i2) && (i2 == i3) && (i3 == i4))
    if x < 0.6
        scatter(x,y,8,'k','filled')
        all_points = all_points+1;
        [~,closest] = min(abs(x-Bayes_noise_values/V1));
        if y <= (TH90(closest)-V0)/V1
           grey_region = grey_region+1;
        end
    end
end
end
end
end
end

% legend

axes_position = [mleft+width_A+gapx1+width_B+gapx3 mbottom+legend_yshift width_legend height_C];
axes('unit','pixel','position',axes_position)
set(gca,'FontName','Helvetica','fontsize',font_size,'linewidth',1)
plot([0 1],[0 .7],'w')
hold on
xlim([0 1])
ylim([0 .7])
axis off

legend_x0 = 0.065; 
scatter(0,.68,200,cell_color,'filled')
text(legend_x0,.685,'measured parameters','fontname','Helvetica','fontsize',font_size-5)
scatter(0,.61,20,'k','filled')
text(legend_x0,.615,'permuted parameters','fontname','Helvetica','fontsize',font_size-5)
scatter(0,.54,200,'s','filled','MarkerFaceColor',grey)
text(legend_x0,.545,'robust performance regime','fontname','Helvetica','fontsize',font_size-5)

set(gcf,'PaperPositionMode','auto','papersize',[37.5 28]);
print(gcf,mfilename,'-dpdf','-r0')
saveas(gcf,[mfilename,'.png']);

all_points
grey_region