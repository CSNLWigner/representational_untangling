load('data_raw.mat');

search_thresholds = [4 2 2 2 3 4 1.2];

figure('units','normalized','outerposition',[0 0 1 1])
for n = 1:size(cells,2)  
    m = sum(strcmp({cells(1:n).type},cells(n).type));
    for trial = 1:cells(n).trials
        eval([sprintf('U_raw = %s%d_trial%d;',cells(n).type,m,trial)]);
        DU = diff(U_raw);
        plot(DU)
        hold on
        line([1 length(DU)],[search_thresholds(n) search_thresholds(n)],'color','r')
        xlim([1 length(DU)]);
        xlabel('time [dt]');
        ylabel('du [mV]');
        title(sprintf('%s %d trial %d',cells(n).type,m,trial));
        hold off
        w = waitforbuttonpress;
        if w == 0
            break
        end
    end
end