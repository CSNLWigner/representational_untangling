ADD = @(id,type,trials) struct('id',id,'type',type,'trials',trials);

cells(1) = ADD('030911_c1_contrast_008_60','simple',6);
cells(2) = ADD('031111_c3_contrast_180_004','simple',5);
cells(3) = ADD('031011_c3_contrast_002_240','simple',5);
cells(4) = ADD('030211_c1_contrast_150','simple',5);
cells(5) = ADD('030211_c2_contrast_150','simple',5);
cells(6) = ADD('031511_c1_contrast_002_330','complex',5);
cells(7) = ADD('032311_c2_contrast_150_002','complex',4);

cd(fileparts(which(mfilename)));
mdatafile = 'data_raw.mat';
save(mdatafile,'cells');

for n = 1:size(cells,2)
    i = sum(strcmp({cells(1:n).type},cells(n).type));
    for trial = 1:cells(n).trials
        load(sprintf('data/data_%s_00%d.mat',cells(n).id,trial));
        varname_data = sprintf('%s%d_trial%d',cells(n).type,i,trial);
        varname_indices = sprintf('%s%d_trial%d_indices',cells(n).type,i,trial);
        eval([sprintf('%s = d(:,2);',varname_data)]);
        indices = [[1;iStart],iEnd]';
        indices = indices(:);
        eval([sprintf('%s = indices(2:end);',varname_indices)]);
        save(mdatafile,varname_data,varname_indices,'-append');
        for contrast = 100:-25:25
            c = contrast/25;
            varname_data = sprintf('%s%d_trial%d_c%d',cells(n).type,i,trial,contrast);
            eval([sprintf('%s = d(iStart(%d):iEnd(%d+1)-1,2);',varname_data,c,c)]);
            save(mdatafile,varname_data,'-append');
        end
        varname_data = sprintf('%s%d_trial%d_unknown1',cells(n).type,i,trial);
        eval([sprintf('%s = d(1:iEnd(1)-1,2);',varname_data)]);
        save(mdatafile,varname_data,'-append');
        for grey = 1:4
            varname_data = sprintf('%s%d_trial%d_grey%d',cells(n).type,i,trial,grey);
            eval([sprintf('%s = d(iEnd(%d):iStart(%d)-1,2);',varname_data,grey,grey)]);
            save(mdatafile,varname_data,'-append');
        end
        varname_data = sprintf('%s%d_trial%d_unknown2',cells(n).type,i,trial);
        eval([sprintf('%s = d(iEnd(5):end,2);',varname_data)]);
        save(mdatafile,varname_data,'-append');
        varname_data = sprintf('%s%d_trial%d_stimulus',cells(n).type,i,trial);
        eval([sprintf('%s = d(iEnd(1):iEnd(5)-1,2);',varname_data)]);
        save(mdatafile,varname_data,'-append');
    end
end

dt = 1000*(d(2,1)-d(1,1));
save(mdatafile,'dt','-append')