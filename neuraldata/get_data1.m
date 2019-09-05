function data = get_data(data_type,cell_index,trial,segment)

varname_trial = sprintf('simple%d_trial%d',cell_index,trial);
varname_indices = sprintf('simple%d_trial%d_indices',cell_index,trial);
load(sprintf('../neuraldata/cell%d_raw.mat',cell_index),varname_trial,varname_indices);

switch segment
    case 'known'
        eval(sprintf('i1 = %s(1);',varname_indices));
        eval(sprintf('i2 = %s(9);',varname_indices));
    case 'unknown1'
        i1 = 1;
        eval(sprintf('i2 = %s(1);',varname_indices));
    case 'unknown2'
        eval(sprintf('i1 = %s(9);',varname_indices));
        eval(sprintf('i2 = length(%s);',varname_trial));
    case 'c100' % contrast 1
        eval(sprintf('i1 = %s(2);',varname_indices));
        eval(sprintf('i2 = %s(3);',varname_indices));
    case 'c75' % contrast 0.75
        eval(sprintf('i1 = %s(4);',varname_indices));
        eval(sprintf('i2 = %s(5);',varname_indices));
    case 'c50' % contrast 0.5
        eval(sprintf('i1 = %s(6);',varname_indices));
        eval(sprintf('i2 = %s(7);',varname_indices));
    case 'c25' % contrast 0.25
        eval(sprintf('i1 = %s(8);',varname_indices));
        eval(sprintf('i2 = %s(9);',varname_indices));  
    case 'grey1' % grey screen
        eval(sprintf('i1 = %s(1);',varname_indices));
        eval(sprintf('i2 = %s(2);',varname_indices));
    case 'grey2' % grey screen
        eval(sprintf('i1 = %s(3);',varname_indices));
        eval(sprintf('i2 = %s(4);',varname_indices));
    case 'grey3' % grey screen
        eval(sprintf('i1 = %s(5);',varname_indices));
        eval(sprintf('i2 = %s(6);',varname_indices));
    case 'grey4' % grey screen
        eval(sprintf('i1 = %s(7);',varname_indices));
        eval(sprintf('i2 = %s(8);',varname_indices));
    otherwise
        i1 = 1;
        eval([sprintf('i2 = length(%s);',varname_trial)]);
end

if strcmp(data_type,'peaks')
    varname_peaks = sprintf('simple%d_trial%d_peaks',cell_index,trial);
    load(sprintf('../neuraldata/cell%d_gen.mat',cell_index),varname_peaks);    
    eval(sprintf('data = %s(%s>=i1 & %s<i2)-i1+1;',varname_peaks,varname_peaks,varname_peaks));
elseif strcmp(data_type,'spindices')
    varname_spindices = sprintf('simple%d_trial%d_spindices',cell_index,trial);
    load(sprintf('../neuraldata/cell%d_gen.mat',cell_index),varname_spindices);    
    eval(sprintf('data = %s(%s>=i1 & %s<i2)-i1+1;',varname_spindices,varname_spindices,varname_spindices));
elseif strcmp(data_type,'thresholds')
    varname_thresholds = sprintf('simple%d_trial%d_thresholds',cell_index,trial);
    varname_peaks = sprintf('%s%d_trial%d_peaks',cell_type,index,trial);
    load(sprintf('../neuraldata/cell%d_gen.mat',cell_index),varname_thresholds,varname_peaks);
    eval(sprintf('data = %s(%s>=i1 & %s<i2);',varname_thresholds,varname_peaks,varname_peaks));
else
    varname_data = sprintf('simple%d_trial%d_%s',cell_index,trial,data_type);
    load(sprintf('../neuraldata/cell%d_gen.mat',cell_index),varname_data);
    eval(sprintf('data = %s(i1:i2-1);',varname_data));
end        

end