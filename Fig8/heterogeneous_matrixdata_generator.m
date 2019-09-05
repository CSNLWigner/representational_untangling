close all
clearvars

cd(fileparts(which(mfilename)));

kappa_mvalues = 1:0.2:1.8;
noise_mvalues = 2:0.5:4;

kappa_msize = length(kappa_mvalues);
noise_msize = length(noise_mvalues);

OPTTH = zeros(kappa_msize,noise_msize);
OPTRELTH = zeros(kappa_msize,noise_msize);

noise_indexshift = 3;

for kappa_index = 1:5

    load(sprintf('kappa1d%d_data.mat',2*(kappa_index-1)),'MAXX','V0','V1')
   
    for i = 1:noise_msize
        i1 = i+noise_indexshift;
        OPTTH(kappa_index,i) = MAXX(i1);
        OPTRELTH(kappa_index,i) = (MAXX(i1)-V0)/V1;
    end
end

save('heterogeneous_matrixdata.mat','kappa_mvalues','noise_mvalues','kappa_msize','noise_msize','V0','V1','OPTTH','OPTRELTH');