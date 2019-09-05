close all
clearvars
cd(fileparts(which(mfilename)));
load('data1_raw.mat')
save('cell1_raw.mat','simple1_trial1','simple1_trial1_indices','simple1_trial2','simple1_trial2_indices','simple1_trial3','simple1_trial3_indices','simple1_trial4','simple1_trial4_indices','simple1_trial5','simple1_trial5_indices','simple1_trial6','simple1_trial6_indices')
save('cell2_raw.mat','simple2_trial1','simple2_trial1_indices','simple2_trial2','simple2_trial2_indices','simple2_trial3','simple2_trial3_indices','simple2_trial4','simple2_trial4_indices','simple2_trial5','simple2_trial5_indices')
save('cell3_raw.mat','simple3_trial1','simple3_trial1_indices','simple3_trial2','simple3_trial2_indices','simple3_trial3','simple3_trial3_indices','simple3_trial4','simple3_trial4_indices','simple3_trial5','simple3_trial5_indices')
save('cell4_raw.mat','simple4_trial1','simple4_trial1_indices','simple4_trial2','simple4_trial2_indices','simple4_trial3','simple4_trial3_indices','simple4_trial4','simple4_trial4_indices','simple4_trial5','simple4_trial5_indices')