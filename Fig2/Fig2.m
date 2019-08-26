% script for creating Fig. 2
% extrenal functions used:
%   cresponse
%   nonlin
%   circular_gabor
%   mixed_rho
%   mixed_p0
%   paramset
close all
clearvars
cd(fileparts(which(mfilename)));
addpath('../functions/');

% predefined colors:
c1 = [0 51 255]/255; % dot blue
c2 = [204 64 32]/255; % dot red
c3 = [1 1 0];
c4 = [0 .85 0];
c5 = [0 1 .8];
c6 = [0 .8 .8];

% generating colors for orientations:
segment_size = [3 7 7 6 6 7];
segment_color = [c1; c2; c3; c4; c5; c6];
segment_color = [segment_color; segment_color(1,:)];
number_of_segments = length(segment_size);
COLORS = [];
for j = 1:number_of_segments
    mix = linspace(1,0,segment_size(j)+1);
    colors = mix'*segment_color(j,:) + (1-mix)'*segment_color(j+1,:);
    COLORS = [COLORS; colors(2:end,:)];
end
COLORS = flipud(COLORS);
COLORS = [COLORS;COLORS];

% line colors:
Bayes_color = [0 .85 0];
dark_green = [0 .6 0];

glambda = 3; gsigma = 2; slambda = glambda;
dgtheta = 30; dgphi = 30;
dstheta = 5; dsphi = 1.25;

% x0, y0, lambda_G, theta_G, phi_G, sigma
gabor1_params = [0,0,glambda,0,0,gsigma];
gabor2_params = [0,0,glambda,dgtheta,dgphi,gsigma];

threshold = .05;
noise_std = .1;

% parameters of stimuli indicated by red and blue dots on panel A
stheta1 = 5; sphi1 = 55;
stheta2 = 15; sphi2 = 40;

s1_params = [slambda,stheta1,sphi1]; % blue
s2_params = [slambda,stheta2,sphi2]; % red

s1_u1 = cresponse(gabor1_params,s1_params);
s1_u2 = cresponse(gabor2_params,s1_params);
s2_u1 = cresponse(gabor1_params,s2_params);
s2_u2 = cresponse(gabor2_params,s2_params);

s1_r1 = nonlin(s1_u1,threshold);
s1_r2 = nonlin(s1_u2,threshold);
s2_r1 = nonlin(s2_u1,threshold);
s2_r2 = nonlin(s2_u2,threshold);

% figure settings:
width = 275;
height = 275;
mtop = 43;
mbottom = 45;
mleft = 48; 
mright = 20;
gapx = 57;
gapy = 75;
figure_width = 3*width+2*gapx+mleft+mright;
figure_height = 2*height+gapy+mbottom+mtop;
set(gcf,'unit','pixel','position',[0 0 figure_width figure_height])
set(gcf,'color','white')

linegap = .01;
axlinewidth = 3;
light_grey = [.8 .8 .8];

rmax = 1;
textsize = 18;
ABC_size = 28;
ABC_x = .075;
ABC_y = 1.06;
linewidth1 = 3;
linewidth2 = 6;
dotsize0 = 10;
dotsize1 = 20;
dotsize2 = 70;
dotsize3 = 150;
N = 75;
% random seeds:
seed1 = 1111;
seed2 = 2000;
seed3 = 3;
seed4 = 4;

% subplot A: MP witzhout noise
axes_position = [mleft mbottom+height+gapy width height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'FontName','Helvetica')

hold on
axis square
axis off
set(gca,'XLim',[-1.05 1.05],'YLim',[-1.05 1.05],'XTick',[],'YTick',[])

% circular arrow for phase
greek_factor = .9;
cx0 = .775;
cy0 = .64;
rr = .1;
tt = linspace(-.65*pi,1.1*pi);
xx = rr*cos(tt)+cx0;
yy = rr*sin(tt)+cy0;
plot(xx,yy,'color','k','linewidth',2)

% arrowhead
ph0 = 1.1*pi;
dph = pi*30/180;
ax0 = cx0+rr*cos(ph0);
ay0 = cy0+rr*sin(ph0);
aw = .055;
tilt = 1.15;
plot([ax0+aw*cos(ph0+dph-tilt*pi/2) ax0 ax0+aw*cos(ph0-dph-tilt*pi/2)],[ay0+aw*sin(ph0+dph-tilt*pi/2) ay0 ay0+aw*sin(ph0-dph-tilt*pi/2)],'color','k','linewidth',2)

text(cx0,cy0+3.6*rr,{'stimulus'},'FontSize',15,'FontName','Helvetica','HorizontalAlign','center');
text(cx0,cy0+2.4*rr,{'phase'},'FontSize',15,'FontName','Helvetica','HorizontalAlign','center');

text(ABC_x,ABC_y,'A','FontName','Helvetica','FontSize',ABC_size,'HorizontalAlign','center','unit','normalized')
text(0,-1.2,'filter response #1','FontName','Helvetica','FontSize',textsize,'HorizontalAlign','center')
texth = text(-1.22,0,'filter response #2','FontName','Helvetica','FontSize',textsize,'HorizontalAlign','center');
set(texth,'Rotation',90);
plot([-1.05 -1.05],[-1.05 1.05],'k','linewidth',axlinewidth)
plot([-1.05 1.05],[-1.05 -1.05],'k','linewidth',axlinewidth)

lx0 = -.65;
ly0 = .2;
lxw = .075;
lyw = .72;

% color curves for orientations
index = 0;
for stheta = 180:-dstheta:dstheta
    index = index+1;
    u1 = [];
    u2 = [];
    for sphi = 0:dsphi:360+dsphi
        wave_params = [slambda,stheta,sphi];
        u1(end+1) = cresponse(gabor1_params,wave_params);
        u2(end+1) = cresponse(gabor2_params,wave_params);
    end
    plot(u1,u2,'color',COLORS(index,:),'linewidth',linewidth1)
    
    % generating colorbar for orientations
    plot([lx0 lx0+lxw],ly0+lyw*[index/36 index/36],'color',COLORS(index+35,:),'linewidth',2.6)
    
    if stheta1 == stheta
        s1_color = COLORS(index,:);
        s1_U1 = u1';
        s1_U2 = u2';
    end
    if stheta2 == stheta
        s2_color = COLORS(index,:);
        s2_U1 = u1';
        s2_U2 = u2';
    end
end

% cardinal directions next to the colorbar:
xgap = .075; lrad = .075;
for ind = 0:9:36
    plot(lx0+lxw+xgap+lrad+lrad*[-1 1]*cos(pi*ind/36),ly0+lyw*[ind/36 ind/36]+lrad*[-1 1]*sin(pi*ind/36),'color',COLORS(ind+35,:),'linewidth',3)
end
text(lx0-.265,ly0+lyw/2,{'stimulus'},'FontSize',15,'FontName','Helvetica','HorizontalAlign','center','rotation',90);
text(lx0-.14,ly0+lyw/2,{'orientation'},'FontSize',15,'FontName','Helvetica','HorizontalAlign','center','rotation',90);

for stheta = dstheta:dstheta:180
    index = index+1;
    u1 = [];
    u2 = [];
    for sphi = 0:dsphi:360+dsphi
        wave_params = [slambda,stheta,sphi];
        u1(end+1) = cresponse(gabor1_params,wave_params);
        u2(end+1) = cresponse(gabor2_params,wave_params);
    end
end
% red and blue dots
scatter(s1_u1,s1_u2,dotsize3,'k','filled')
scatter(s2_u1,s2_u2,dotsize3,'k','filled')
scatter(s1_u1,s1_u2,dotsize2,s1_color,'filled')
scatter(s2_u1,s2_u2,dotsize2,s2_color,'filled')

gwidth = .1*width;
axes_position = [mleft+width-1.5*gwidth mbottom+height+gapy-1.25*gwidth gwidth gwidth];
subplot_axes = axes('unit','pixel','position',axes_position);
hold on
axis square
axis off

% making a Gabor patch for filter 1
NN = 1000;
G1 = circular_gabor([5],gabor1_params,NN);
maxg = max(max(G1));
for i = 1:NN
    for j = 1:NN
        if (i-NN/2)^2+(j-NN/2)^2 > (NN/2)^2
            G1(i,j) = maxg;
        end
    end 
end  
imagesc(G1)
colormap(.25+.75*gray(256));
cx0 = NN/2;
cy0 = NN/2;
rr = NN/2-1;
tt = linspace(0,2*pi);
xx = cx0+rr*cos(tt);
yy = cy0+rr*sin(tt);
hold on
plot(xx,yy,'color','k','linewidth',1)

axes_position = [mleft-1.25*gwidth mbottom+height+gapy+width-1.5*gwidth gwidth gwidth];
subplot_axes = axes('unit','pixel','position',axes_position);
hold on
axis square
axis off

% making a Gabor patch for filter 2
NN = 1000;
G1 = circular_gabor([5],gabor2_params,NN);
maxg = max(max(G1));
for i = 1:NN
    for j = 1:NN
        if (i-NN/2)^2+(j-NN/2)^2 > (NN/2)^2
            G1(i,j) = maxg;
        end
    end 
end  
imagesc(G1)
colormap(.25+.75*gray(256));
cx0 = NN/2;
cy0 = NN/2;
rr = NN/2-1;
tt = linspace(0,2*pi);
xx = cx0+rr*cos(tt);
yy = cy0+rr*sin(tt);
hold on
plot(xx,yy,'color','k','linewidth',1)

% subplot D: FR without noise
axes_position = [mleft mbottom width height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'FontName','Helvetica')
hold on
axis square
axis off
set(gca,'XLim',[-linegap rmax],'YLim',[-linegap rmax],'XTick',[],'YTick',[])
text(ABC_x,ABC_y,'D','FontSize',ABC_size,'HorizontalAlign','center','unit','normalized')
text(.5,-.075,'truncated filter response #1','FontName','Helvetica','FontSize',textsize,'HorizontalAlign','center')
texth = text(-.09,.5,'truncated filter response #2','FontName','Helvetica','FontSize',textsize,'HorizontalAlign','center');
set(texth,'Rotation',90);
plot([-linegap -linegap],[-linegap 1],'k','linewidth',axlinewidth)
plot([-linegap 1],[-linegap -linegap],'k','linewidth',axlinewidth)
index = 0;
for stheta = 180:-dstheta:dstheta
    index = index+1;
    r1 = [];
    r2 = [];
    for sphi = 0:dsphi:360+dsphi
        wave_params = [slambda,stheta,sphi];
        r1(end+1) = nonlin(cresponse(gabor1_params,wave_params),threshold);
        r2(end+1) = nonlin(cresponse(gabor2_params,wave_params),threshold);
    end
    plot(r1,r2,'color',COLORS(index,:),'linewidth',linewidth1) 
end
for stheta = dstheta:dstheta:180
    index = index+1;
    r1 = [];
    r2 = [];
    for sphi = 0:dsphi:360+dsphi
        wave_params = [slambda,stheta,sphi];
        r1(end+1) = nonlin(cresponse(gabor1_params,wave_params),threshold);
        r2(end+1) = nonlin(cresponse(gabor2_params,wave_params),threshold);
    end
end
scatter(s1_r1,s1_r2,dotsize3,'k','filled')
scatter(s2_r1,s2_r2,dotsize3,'k','filled')
scatter(s1_r1,s1_r2,dotsize2,s1_color,'filled')
scatter(s2_r1,s2_r2,dotsize2,s2_color,'filled')

% subplot B: MP with noise (no nuisance)
axes_position = [mleft+width+gapx mbottom+height+gapy width height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'FontName','Helvetica')
hold on
axis square
axis off
set(gca,'XLim',[-1.05 1.05],'YLim',[-1.05 1.05],'XTick',[],'YTick',[])
text(ABC_x,ABC_y,'B','FontSize',ABC_size,'HorizontalAlign','center','unit','normalized')
text(0,1.18,'fixed phase','FontName','Helvetica','FontSize',textsize,'HorizontalAlignment','center')
text(0,-1.2,'membrane potential response #1','FontName','Helvetica','FontSize',textsize,'HorizontalAlign','center')
texth = text(-1.22,0,'membrane potential response #2','FontName','Helvetica','FontSize',textsize,'HorizontalAlign','center');
set(texth,'Rotation',90);
plot([-1.05 -1.05],[-1.05 1.05],'k','linewidth',axlinewidth)
plot([-1.05 1.05],[-1.05 -1.05],'k','linewidth',axlinewidth)
for stheta = dstheta:dstheta:180
    u1 = [];
    u2 = [];
    for sphi = 0:dsphi:360+dsphi
        wave_params = [slambda,stheta,sphi];
        u1(end+1) = cresponse(gabor1_params,wave_params);
        u2(end+1) = cresponse(gabor2_params,wave_params);
    end
    plot(u1,u2,'color',light_grey,'linewidth',linewidth1) 
end
rng(seed1)
s1_u1_noisy = s1_u1+noise_std*randn(N,1);
s1_u2_noisy = s1_u2+noise_std*randn(N,1);
rng(seed2)
s2_u1_noisy = s2_u1+noise_std*randn(N,1);
s2_u2_noisy = s2_u2+noise_std*randn(N,1);
scatter(s1_u1_noisy,s1_u2_noisy,dotsize1,s1_color,'filled')
scatter(s2_u1_noisy,s2_u2_noisy,dotsize1,s2_color,'filled')
scatter(s1_u1,s1_u2,dotsize3,'k','filled')
scatter(s2_u1,s2_u2,dotsize3,'k','filled')
scatter(s1_u1,s1_u2,dotsize2,light_grey,'filled')
scatter(s2_u1,s2_u2,dotsize2,light_grey,'filled')

MP_threshold = -5;

NN = 100;
xspace = linspace(-1,1,NN);
yspace = linspace(-1,1,NN);
[XCOORD,YCOORD] = meshgrid(xspace,yspace);
Z = zeros(NN,NN);

GBANK = [gabor1_params; gabor2_params];
SBANK0 = [s1_params; s2_params];

MU = 0+1*cresponse(GBANK,SBANK0); % mu_n(theta_k,phi_j) = MU((k-1)*J+j,n)

PK_OPT = zeros(2,NN^2);

ii = 0;
for ix = 1:NN
    for iy = 1:NN
        ii = ii+1;
        r1 = nonlin(xspace(ix),MP_threshold,1);
        r2 = nonlin(yspace(iy),MP_threshold,1);
        R = [r1 r2];
        for kk = 1:2
            product = 1;
            for nn = 1:2
                if R(nn) == 0
                    product = product*mixed_p0(MU(kk,nn),noise_std,MP_threshold);
                else
                    product = product*mixed_rho(R(nn),MU(kk,nn),noise_std,MP_threshold,1,1);
                end
            end
            PK_OPT(kk,ii) = PK_OPT(kk,ii) + product;
        end
        PK_OPT(:,ii) = PK_OPT(:,ii)/sum(PK_OPT(:,ii));
        Z(ix,iy) = PK_OPT(1,ii);
    end
end
contour(YCOORD,XCOORD,Z,[.5 .5],'linewidth',3,'color',Bayes_color)

% subplot E: FR with noise (no nuisance)
axes_position = [mleft+width+gapx mbottom width height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'FontName','Helvetica')
hold on
axis square
axis off
set(gca,'XLim',[-linegap rmax],'YLim',[-linegap rmax],'XTick',[],'YTick',[])
text(ABC_x,ABC_y,'E','FontSize',ABC_size,'HorizontalAlign','center','unit','normalized')
text(.5,-.075,'firing rate response #1','FontName','Helvetica','FontSize',textsize,'HorizontalAlign','center')
texth = text(-.09,.5,'firing rate response #2','FontName','Helvetica','FontSize',textsize,'HorizontalAlign','center');
set(texth,'Rotation',90);
plot([-linegap -linegap],[-linegap 1],'k','linewidth',axlinewidth)
plot([-linegap 1],[-linegap -linegap],'k','linewidth',axlinewidth)
for stheta = dstheta:dstheta:180
    r1 = [];
    r2 = [];
    for sphi = 0:dsphi:360+dsphi
        wave_params = [slambda,stheta,sphi];
        r1(end+1) = nonlin(cresponse(gabor1_params,wave_params),threshold);
        r2(end+1) = nonlin(cresponse(gabor2_params,wave_params),threshold);
    end
    plot(r1,r2,'color',light_grey,'linewidth',linewidth1) 
end
scatter(nonlin(s1_u1_noisy,threshold),nonlin(s1_u2_noisy,threshold),dotsize1,s1_color,'filled')
scatter(nonlin(s2_u1_noisy,threshold),nonlin(s2_u2_noisy,threshold),dotsize1,s2_color,'filled')
scatter(s1_r1,s1_r2,dotsize3,'k','filled')
scatter(s2_r1,s2_r2,dotsize3,'k','filled')
scatter(s1_r1,s1_r2,dotsize2,light_grey,'filled')
scatter(s2_r1,s2_r2,dotsize2,light_grey,'filled')

ii = 0;
for ix = 1:NN
    for iy = 1:NN
        ii = ii+1;
        r1 = nonlin(xspace(ix),threshold,1);
        r2 = nonlin(yspace(iy),threshold,1);
        R = [r1 r2];    
        for kk = 1:2
            product = 1;
            for nn = 1:2
                if R(nn) == 0
                    product = product*mixed_p0(MU(kk,nn),noise_std,threshold);
                else
                    product = product*mixed_rho(R(nn),MU(kk,nn),noise_std,threshold,1,1);
                end
            end
            PK_OPT(kk,ii) = PK_OPT(kk,ii) + product;    
        end
        PK_OPT(:,ii) = PK_OPT(:,ii)/sum(PK_OPT(:,ii));
        Z(ix,iy) = PK_OPT(1,ii);
    end
end
contour(YCOORD,XCOORD,Z,[.5 .5],'linewidth',3,'color',Bayes_color)

% subplot C: MP with noise (with nuisance)
axes_position = [mleft+2*width+2*gapx mbottom+height+gapy width height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'FontName','Helvetica')
hold on
axis square
axis off
set(gca,'XLim',[-1.05 1.05],'YLim',[-1.05 1.05],'XTick',[],'YTick',[])
text(ABC_x,ABC_y,'C','FontSize',ABC_size,'HorizontalAlign','center','unit','normalized')
text(0,1.18,'variable phase','FontName','Helvetica','FontSize',textsize,'HorizontalAlignment','center')
plot([-1.05 -1.05],[-1.05 1.05],'k','linewidth',axlinewidth)
plot([-1.05 1.05],[-1.05 -1.05],'k','linewidth',axlinewidth)
text(0,-1.2,'membrane potential response #1','FontName','Helvetica','FontSize',textsize,'HorizontalAlign','center')
texth = text(-1.22,0,'membrane potential response #2','FontName','Helvetica','FontSize',textsize,'HorizontalAlign','center');
set(texth,'Rotation',90);
for stheta = dstheta:dstheta:180
    u1 = [];
    u2 = [];
    for sphi = 0:dsphi:360+dsphi
        wave_params = [slambda,stheta,sphi];
        u1(end+1) = cresponse(gabor1_params,wave_params);
        u2(end+1) = cresponse(gabor2_params,wave_params);
    end
    plot(u1,u2,'color',light_grey,'linewidth',linewidth1)
end
for stheta = dstheta:dstheta:180
    u1 = [];
    u2 = [];
    for sphi = 0:dsphi:360+dsphi
        wave_params = [slambda,stheta,sphi];
        u1(end+1) = cresponse(gabor1_params,wave_params);
        u2(end+1) = cresponse(gabor2_params,wave_params);
    end
    if (stheta == stheta1) || (stheta == stheta2)
        plot(u1,u2,'color','k','linewidth',linewidth2)
        plot(u1,u2,'color',light_grey,'linewidth',linewidth1)
    end
end
M = 1;
s1_U1 = repmat(s1_U1,M,1)+noise_std*randn(M*length(s1_U1),1);
s1_U2 = repmat(s1_U2,M,1)+noise_std*randn(M*length(s1_U2),1);
s2_U1 = repmat(s2_U1,M,1)+noise_std*randn(M*length(s2_U1),1);
s2_U2 = repmat(s2_U2,M,1)+noise_std*randn(M*length(s2_U2),1);
scatter(s1_U1,s1_U2,dotsize1,s1_color,'filled')
scatter(s2_U1,s2_U2,dotsize1,s2_color,'filled')

%

J = 50;
s0lambda = {'const',slambda};
s1theta = {'const',stheta1};
s2theta = {'const',stheta2};
s0phi = {'grid',0,360,J};
SBANK1 = paramset(s0lambda,s1theta,s0phi);
SBANK2 = paramset(s0lambda,s2theta,s0phi);
SBANK0 = [SBANK1; SBANK2];
Y0 = [ones(J,1); 2*ones(J,1)];

MU = 0+1*cresponse(GBANK,SBANK0); % mu_n(theta_k,phi_j) = MU((k-1)*J+j,n)

PK_OPT = zeros(2,NN^2);

% optimal decoding
ii = 0;
for ix = 1:NN
    for iy = 1:NN
        ii = ii+1;
        r1 = nonlin(xspace(ix),MP_threshold,1);
        r2 = nonlin(yspace(iy),MP_threshold,1);
        R = [r1 r2];
        for kk = 1:2
            for jj = 1:J
                product = 1;
                for nn = 1:2
                    if R(nn) == 0
                        product = product*mixed_p0(MU((kk-1)*J+jj,nn),noise_std,MP_threshold);
                    else
                        product = product*mixed_rho(R(nn),MU((kk-1)*J+jj,nn),noise_std,MP_threshold,1,1);
                    end
                end
                PK_OPT(kk,ii) = PK_OPT(kk,ii) + product;
            end 
        end
        PK_OPT(:,ii) = PK_OPT(:,ii)/sum(PK_OPT(:,ii));
        Z(ix,iy) = PK_OPT(1,ii);
    end
end
contour(YCOORD,XCOORD,Z,[.5 .5],'linewidth',3,'color',Bayes_color)

% subplot F: FR with noise (with nuisance)
axes_position = [mleft+2*width+2*gapx mbottom width height];
subplot_axes = axes('unit','pixel','position',axes_position);
set(gca,'FontName','Helvetica')
hold on
axis square
axis off
set(gca,'XLim',[-linegap rmax],'YLim',[-linegap rmax],'XTick',[],'YTick',[])
text(ABC_x,ABC_y,'F','FontSize',ABC_size,'HorizontalAlign','center','unit','normalized')
text(.5,-.075,'firing rate response #1','FontName','Helvetica','FontSize',textsize,'HorizontalAlign','center')
texth = text(-.09,.5,'firing rate response #2','FontName','Helvetica','FontSize',textsize,'HorizontalAlign','center');
set(texth,'Rotation',90);
plot([-linegap -linegap],[-linegap 1],'k','linewidth',axlinewidth)
plot([-linegap 1],[-linegap -linegap],'k','linewidth',axlinewidth)
for stheta = dstheta:dstheta:180
    r1 = [];
    r2 = [];
    for sphi = 0:dsphi:360+dsphi
        wave_params = [slambda,stheta,sphi];
        r1(end+1) = nonlin(cresponse(gabor1_params,wave_params),threshold);
        r2(end+1) = nonlin(cresponse(gabor2_params,wave_params),threshold);
    end
    plot(r1,r2,'color',light_grey,'linewidth',linewidth1) 
end
for stheta = dstheta:dstheta:180
    r1 = [];
    r2 = [];
    for sphi = 0:dsphi:360+dsphi
        wave_params = [slambda,stheta,sphi];
        r1(end+1) = nonlin(cresponse(gabor1_params,wave_params),threshold);
        r2(end+1) = nonlin(cresponse(gabor2_params,wave_params),threshold);
    end
    if (stheta == stheta1) || (stheta == stheta2)
        plot(r1,r2,'color','k','linewidth',linewidth2)
        plot(r1,r2,'color',light_grey,'linewidth',linewidth1)
    end
end

scatter(nonlin(s1_U1,threshold),nonlin(s1_U2,threshold),20,s1_color,'filled')
scatter(nonlin(s2_U1,threshold),nonlin(s2_U2,threshold),20,s2_color,'filled')

% optimal decoding
ii = 0;
for ix = 1:NN
    for iy = 1:NN
        ii = ii+1;
        r1 = nonlin(xspace(ix),threshold,1);
        r2 = nonlin(yspace(iy),threshold,1);
        R = [r1 r2];
        for kk = 1:2
            for jj = 1:J
                product = 1;
                for nn = 1:2
                    if R(nn) == 0
                        product = product*mixed_p0(MU((kk-1)*J+jj,nn),noise_std,threshold);
                    else
                        product = product*mixed_rho(R(nn),MU((kk-1)*J+jj,nn),noise_std,threshold,1,1);
                    end
                end
                PK_OPT(kk,ii) = PK_OPT(kk,ii) + product;
            end 
        end
        PK_OPT(:,ii) = PK_OPT(:,ii)/sum(PK_OPT(:,ii));
        Z(ix,iy) = PK_OPT(1,ii);
    end
end
plot([0 1],[.03 .46],'linewidth',3,'color',dark_green)
contour(YCOORD,XCOORD,Z,[.5 .5],'linewidth',3,'color',Bayes_color)

set(gcf,'PaperPositionMode','auto','papersize',[36 26]);
saveas(gcf,sprintf('%s.pdf',mfilename));