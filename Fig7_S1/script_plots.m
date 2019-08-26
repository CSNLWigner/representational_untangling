clear all
close all

load('script_2d_privatenoise_ref_data.mat','thresholds','FC_LIN','PFC_LIN');
thresholds1 = thresholds;
FC_REF = FC_LIN;
PFC_REF = PFC_LIN;

load('script_2d_ILC_data.mat','FC_LIN','PFC_LIN');
FC_ILC = FC_LIN;
PFC_ILC = PFC_LIN;

load('script_2d_Bayes_ILC_data.mat','thresholds','FC_OPT','PFC_OPT');

figure('units','normalized','outerposition',[0 0 1 1],'color','w')

subplot(1,2,1)
hold on

p1 = plot(thresholds1,FC_REF,'color','b','linewidth',2);
p2 = plot(thresholds1,FC_ILC,'color','k','linewidth',2);
p3 = plot(thresholds,FC_OPT,'color','r','linewidth',2);

legend([p1,p2,p3],{'ref','linear ILC','Bayes'})

axis([thresholds(1) thresholds(end) 0 1])
xlabel('threshold [mV]');
ylabel('fraction correct');
title('FC')

subplot(1,2,2)
hold on

p1 = plot(thresholds1,PFC_REF,'color','b','linewidth',2);
p2 = plot(thresholds1,PFC_ILC,'color','k','linewidth',2);
p3 = plot(thresholds,PFC_OPT,'color','r','linewidth',2);

legend([p1,p2,p3],{'ref','linear ILC','Bayes'})

axis([thresholds(1) thresholds(end) 0 1])
xlabel('threshold [mV]');
ylabel('probabilistic fraction correct');
title('PFC')

saveas(gcf,'plots.png');