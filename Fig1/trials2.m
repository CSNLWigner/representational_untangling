% script for generating part of Fig. 1
% MP and FR illustrations for cell N (trials and their mean)
close all
clearvars
set(gcf,'color','white')
seeds = [1030,1008,1015,1000];
K = length(seeds); % number_of_trials
lws = 1; % square box linewidth
lwb = 4; % black linewidth
lwg = 3.5; % grey linewidth
lambda = 1; % wavelength (unit of time period)
cycles = 1.2; % number of cycles
L = 250; % length of time vector
time = linspace(0,cycles*lambda,L);
AC = .5; % amplitude
maxy = 1.25;
phi = -30; % phase
DC = .2; 
Y = DC+AC*sin(2*pi*time/lambda+phi);
noise = .5;
trials = repmat(Y,K,1);
kernel_sigma = .025;
for k = K:-1:1
    rng(seeds(k));
    trials(k,:) = trials(k,:)+noise*randn(1,L);
    for i = 1:L
        KERNEL = exp(-(time-time(i)).^2/kernel_sigma^2/2);
        KERNEL = KERNEL/sum(KERNEL);
        trials(k,i) = sum(KERNEL.*trials(k,:));
    end
    if k == 1
        lcolor = [0 0 0];
        lwidth = lwb; 
    else
        lcolor = [.75 .75 .75];
        lwidth = lwg;
    end  
    plot(time,trials(k,:),'color',lcolor,'linewidth',lwidth)
    hold on
end
ylim([-maxy,maxy])
axis square
set(gca,'xtick',[],'ytick',[],'linewidth',lws,'layer','top')
cd(fileparts(which(mfilename)));
saveas(gcf,'trials2_MP.pdf');
close all
set(gcf,'color','white')
threshold = DC+0.25;
rates = trials-threshold;
rates(rates<0) = 0;
epsilon = .025;
factor = 1.25;
for k = K:-1:1
    if k == 1
        lcolor = [0 0 0];
        lwidth = lwb; 
    else
        lcolor = [.75 .75 .75];
        lwidth = lwg;
    end  
    plot(time,epsilon+factor*rates(k,:),'color',lcolor,'linewidth',lwidth)
    hold on
end
ylim([0,maxy])
axis square
set(gca,'xtick',[],'ytick',[],'linewidth',lws,'layer','top')
saveas(gcf,'trials2_FR.pdf');