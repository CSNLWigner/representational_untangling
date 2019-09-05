function data = get_data1(data_type,cell_index,trial,segment)

varname_trial = sprintf('simple%d_trial%d',cell_index,trial);
varname_indices = sprintf('simple%d_trial%d_indices',cell_index,trial);
load(sprintf('../neuraldata/cell%d_raw.mat',cell_index),varname_trial,varname_indices);

switch segment
    case 'known'
        eval([sprintf('i1 = %s(1);',varname_indices)]);
        eval([sprintf('i2 = %s(9);',varname_indices)]);
    case 'unknown1'
        i1 = 1;
        eval([sprintf('i2 = %s(1);',varname_indices)]);
    case 'unknown2'
        eval([sprintf('i1 = %s(9);',varname_indices)]);
        eval([sprintf('i2 = length(%s);',varname_trial)]);
    case 'c100'
        eval([sprintf('i1 = %s(2);',varname_indices)]);
        eval([sprintf('i2 = %s(3);',varname_indices)]);
    case 'c75'
        eval([sprintf('i1 = %s(4);',varname_indices)]);
        eval([sprintf('i2 = %s(5);',varname_indices)]);
    case 'c50'
        eval([sprintf('i1 = %s(6);',varname_indices)]);
        eval([sprintf('i2 = %s(7);',varname_indices)]);
    case 'c25'
        eval([sprintf('i1 = %s(8);',varname_indices)]);
        eval([sprintf('i2 = %s(9);',varname_indices)]);  
    case 'grey1'
        eval([sprintf('i1 = %s(1);',varname_indices)]);
        eval([sprintf('i2 = %s(2);',varname_indices)]);
    case 'grey2'
        eval([sprintf('i1 = %s(3);',varname_indices)]);
        eval([sprintf('i2 = %s(4);',varname_indices)]);
    case 'grey3'
        eval([sprintf('i1 = %s(5);',varname_indices)]);
        eval([sprintf('i2 = %s(6);',varname_indices)]);
    case 'grey4'
        eval([sprintf('i1 = %s(7);',varname_indices)]);
        eval([sprintf('i2 = %s(8);',varname_indices)]);
    otherwise
        i1 = 1;
        eval([sprintf('i2 = length(%s);',varname_trial)]);
end

if strcmp(data_type,'peaks')
    varname_peaks = sprintf('%s%d_trial%d_peaks',cell_type,index,trial);
    load('data1_gen.mat',varname_peaks);    
    eval([sprintf('data = %s(%s>=i1 & %s<i2)-i1+1;',varname_peaks,varname_peaks,varname_peaks)]);
elseif strcmp(data_type,'spindices')
    varname_spindices = sprintf('%s%d_trial%d_spindices',cell_type,index,trial);
    load('data1_gen.mat',varname_spindices);    
    eval([sprintf('data = %s(%s>=i1 & %s<i2)-i1+1;',varname_spindices,varname_spindices,varname_spindices)]);
elseif strcmp(data_type,'thresholds')
    varname_thresholds = sprintf('%s%d_trial%d_thresholds',cell_type,index,trial);
    varname_peaks = sprintf('%s%d_trial%d_peaks',cell_type,index,trial);
    load('data1_gen.mat',varname_thresholds,varname_peaks);
    eval([sprintf('data = %s(%s>=i1 & %s<i2);',varname_thresholds,varname_peaks,varname_peaks)]);
else
    varname_data = sprintf('%s%d_trial%d_%s',cell_type,index,trial,data_type);
    load('data1_gen.mat',varname_data);
    eval([sprintf('data = %s(i1:i2-1);',varname_data)]);
end        

end