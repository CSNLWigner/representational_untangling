clear all
close all

nonlinth = [-44.8729658328067;-38.6225420909711;-44.8735910204191;-35.0601801865569];
nonlinth_c100 = [-44.7478874740286;-38.4811015151969;-44.7471553930048;-36.3470263609646];

cell_colors = [1 .3 .3; .4 .4 1; 1 .5 0; .2 .7 .5];

% rule 1: homogeneous quality
% rule 2: minumum duration = 4 cycle

% cell 1

segment(1).cell_index = 1;
segment(1).trial_index = 2;
segment(1).start = 0;
segment(1).duration = 6;
segment(1).flag = '';

segment(2).cell_index = 1;
segment(2).trial_index = 3;
segment(2).start = 0;
segment(2).duration = 6;
segment(2).flag = '';

segment(3).cell_index = 1;
segment(3).trial_index = 4;
segment(3).start = 0;
segment(3).duration = 6;
segment(3).flag = 'odd one out';

segment(4).cell_index = 1;
segment(4).trial_index = 5;
segment(4).start = 0;
segment(4).duration = 5;
segment(4).flag = '';

% cell 2

segment(5).cell_index = 2;
segment(5).trial_index = 1;
segment(5).start = 0;
segment(5).duration = 4;
segment(5).flag = '';

segment(6).cell_index = 2;
segment(6).trial_index = 3;
segment(6).start = 0;
segment(6).duration = 6;
segment(6).flag = '';

% cell 3

segment(7).cell_index = 3;
segment(7).trial_index = 3;
segment(7).start = 1;
segment(7).duration = 4;
segment(7).flag = '';

segment(8).cell_index = 3;
segment(8).trial_index = 5;
segment(8).start = 1;
segment(8).duration = 4;
segment(8).flag = 'bad quality';

% cell 4

segment(9).cell_index = 4;
segment(9).trial_index = 1;
segment(9).start = 0.5;
segment(9).duration = 4;
segment(9).flag = 'bad quality';

cd(fileparts(which(mfilename)));
load('data1_raw.mat','dt');
dt = dt/1000;
msec = 1/1000;
tbin = 20*msec; 
nbin = round(tbin/dt);

trials = [6 5 5 5];
TC = 0.5; % sec
for i = 1:length(segment)
    
    ugen_c100 = get_data1('genpot','simple',segment(i).cell_index,segment(i).trial_index,'c100');
    N = length(ugen_c100);
    
    n1 = round(N*segment(i).start/6)+1;
    n2 = round(N*(segment(i).start+segment(i).duration)/6);
    u_fit = ugen_c100(n1:n2)';
    t_fit = linspace(0,segment(i).duration*TC,length(u_fit));
    
    minu = min(u_fit);
    maxu = max(u_fit);
    y_DC = mean(u_fit);
    y_AC = (maxu-minu)/2;
    
    F1 = @(p,t) p(1)+p(2)*sin(2*pi*t/TC+p(3));
    hint1 = [y_DC,y_AC,0];
    fitoptions = optimoptions('lsqcurvefit','display','off');
    params1 = lsqcurvefit(F1,hint1,t_fit,u_fit,[],[],fitoptions);
    Y1 = F1(params1,t_fit);
    
    F2 = @(p,t) params1(1)+params1(2)*(1-power(2,1-p(1))*power(1+cos(2*pi*t/TC+params1(3)+pi/2),p(1)));
    hint2 = [1.5];
    params2 = lsqcurvefit(F2,hint2,t_fit,u_fit,[],[],fitoptions);
    Y2 = F2(params2,t_fit);

    F3 = @(p,t) p(1)+p(2)*(1-power(2,1-params2(1))*power(1+cos(2*pi*t/TC+params1(3)+pi/2),params2(1)));
    hint3 = [params1(1) params1(2)];
    params3 = lsqcurvefit(F3,hint3,t_fit,u_fit,[],[],fitoptions);
    Y3 = F3(params3,t_fit);

    segment(i).hint_DC = hint1(1);
    segment(i).sinfit_DC = params1(1); 
    segment(i).sinfit_AC = abs(params1(2));
    segment(i).sinfit3_DC = params3(1); 
    segment(i).sinfit3_AC = abs(params3(2));
    
    residum_Y3 = u_fit-Y3;
    lower_residum = residum_Y3(Y3 < params3(1));
    residum_20ms = [];
    for j = 1:floor(length(lower_residum)/nbin)
        range = (j-1)*nbin+1:j*nbin;
        res_mean = mean(lower_residum(range));
        residum_20ms = [residum_20ms,res_mean];
    end

    segment(i).lowernoise_20ms = std(residum_20ms);

%     hold off
%     plot(t_fit,u_fit,'color',[.5 .5 .5 ])
%     hold on
%     plot(t_fit,Y3,'r','linewidth',5)
%     w = waitforbuttonpress;
%     if w == 0
%         break
%     end
    
    cycles_20ms = [];
    for k = 1:segment(i).duration
        
        res = mod(n2-n1+1,segment(i).duration);
        nn = n2-n1+1-res;
        nk = nn/segment(i).duration;
        uk = u_fit((k-1)*nk+1:k*nk);
        
        uk_20ms = [];
        for j = 1:floor(length(uk)/nbin)
            range = (j-1)*nbin+1:j*nbin;
            bin_mean = mean(uk(range));
            uk_20ms = [uk_20ms,bin_mean];
        end
        
        cycles_20ms = [cycles_20ms; uk_20ms];
    end
    
    umean_20ms = mean(cycles_20ms);
    cycles_res20ms = cycles_20ms-repmat(umean_20ms,segment(i).duration,1);
    cycles_res20ms = cycles_res20ms(:,umean_20ms < y_DC);
    
    segment(i).lowerstd_20ms = std(cycles_res20ms(:));
    
    ugen_grey4 = get_data1('genpot','simple',segment(i).cell_index,segment(i).trial_index,'grey4');
    ugen50_grey4 = ugen_grey4(round(length(ugen_grey4)/2):end);
    
    binned_grey4 = [];
    for j = 1:floor(length(ugen_grey4)/nbin)
        range = (j-1)*nbin+1:j*nbin;
        bin_mean = mean(ugen_grey4(range));
        binned_grey4 = [binned_grey4,bin_mean];
    end
    binned50_grey4 = [];
    for j = 1:floor(length(ugen50_grey4)/nbin)
        range = (j-1)*nbin+1:j*nbin;
        bin_mean = mean(ugen50_grey4(range));
        binned50_grey4 = [binned50_grey4,bin_mean];
    end
    
    segment(i).grey_std20ms = std(binned_grey4);
    segment(i).grey50_std20ms = std(binned50_grey4);
    
    segment(i).nonlinth = nonlinth(segment(i).cell_index);
    segment(i).nonlinth_c100 = nonlinth_c100(segment(i).cell_index);
end

save('c100_segments.mat','segment')

subplot(2,2,1)
hold on
for i = 1:length(segment)
    scatter(segment(i).lowernoise_20ms,segment(i).grey_std20ms,100,cell_colors(segment(i).cell_index,:),'filled')
end
xlabel('lower noise /20ms/')
ylabel('grey std /20ms/')
params = polyfit(cell2mat({segment.lowernoise_20ms}),cell2mat({segment.grey_std20ms}),1);
title({sprintf('intercept = %g',params(2)) sprintf('slope = %g',params(1))})

subplot(2,2,2)
hold on
for i = 1:length(segment)
    scatter(segment(i).lowerstd_20ms,segment(i).grey_std20ms,100,cell_colors(segment(i).cell_index,:),'filled')
end
xlabel('lower std /20ms/')
ylabel('grey std /20ms/')
params = polyfit(cell2mat({segment.lowerstd_20ms}),cell2mat({segment.grey_std20ms}),1);
title({sprintf('intercept = %g',params(2)) sprintf('slope = %g',params(1))})

subplot(2,2,3)
hold on
for i = 1:length(segment)
    scatter(segment(i).lowernoise_20ms,segment(i).grey50_std20ms,100,cell_colors(segment(i).cell_index,:),'filled')
end
xlabel('lower noise /20ms/')
ylabel('grey50 std /20ms/')
params = polyfit(cell2mat({segment.lowernoise_20ms}),cell2mat({segment.grey50_std20ms}),1);
title({sprintf('intercept = %g',params(2)) sprintf('slope = %g',params(1))})

subplot(2,2,4)
hold on
for i = 1:length(segment)
    scatter(segment(i).lowerstd_20ms,segment(i).grey50_std20ms,100,cell_colors(segment(i).cell_index,:),'filled')
end
xlabel('lower std /20ms/')
ylabel('grey50 std /20ms/')
params = polyfit(cell2mat({segment.lowerstd_20ms}),cell2mat({segment.grey50_std20ms}),1);
title({sprintf('intercept = %g',params(2)) sprintf('slope = %g',params(1))})