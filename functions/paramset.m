function bank = paramset(varargin)
% Parameter set generator (for making a stimulus bank)
% Input variable syntax: {param_type,param1,param2,...}
% Accepted types: const, grid, rand, unirand, perm, shift, mix
% Example: paramset_generator({'const',1},{'grid',0,360,100});
totaldim = 1;
for k = 1:nargin;
    if findstr(varargin{k}{1},'grid')
        totaldim = totaldim*varargin{k}{4};
    end
end
bank = zeros(totaldim,nargin);
subdim = totaldim;
subrep = 1;
for k = 1:nargin;
    if strcmp(varargin{k}{1},'const')
        bank(:,k) = varargin{k}{2};
    end
    if strcmp(varargin{k}{1},'grid')
        delta = (varargin{k}{3}-varargin{k}{2})/varargin{k}{4};
        grid = varargin{k}{2}+delta/2:delta:varargin{k}{3}-delta/2;
        subdim = subdim/varargin{k}{4};
        repgrid = repmat(grid,subdim,1);
        bank(:,k) = repmat(repgrid(:),subrep,1);
        subrep = subrep*varargin{k}{4};  
    end
    if strcmp(varargin{k}{1},'rand')
        rng(varargin{k}{end});
        bank(:,k) = varargin{k}{2}+(varargin{k}{3}-varargin{k}{2})*rand(totaldim,1);
    end
    if strcmp(varargin{k}{1},'unirand')
        rng(varargin{k}{end});
        delta = (varargin{k}{3}-varargin{k}{2})/totaldim;
        left = varargin{k}{2}:delta:varargin{k}{3}-delta;
        bank(:,k) = left+delta*rand(totaldim,1);
    end
    if strcmp(varargin{k}{1},'perm')
        rng(varargin{k}{end});
        delta = (varargin{k}{3}-varargin{k}{2})/totaldim;
        grid = varargin{k}{2}+delta/2:delta:varargin{k}{3}-delta/2;
        bank(:,k) = grid(randperm(totaldim));    
    end
    if strcmp(varargin{k}{1},'shift')
        delta = (varargin{k}{3}-varargin{k}{2})/totaldim;
        grid = varargin{k}{2}+delta/2:delta:varargin{k}{3}-delta/2;
        order = [1:totaldim,1:totaldim];
        bank(:,k) = grid(order(1+varargin{k}{4}:totaldim+varargin{k}{4}));
    end
    if strcmp(varargin{k}{1},'mix')
        delta = (varargin{k}{3}-varargin{k}{2})/totaldim;
        grid = varargin{k}{2}+delta/2:delta:varargin{k}{3}-delta/2;
        order = 1:totaldim;
        shift = order+ceil(totaldim/2);
        mix = [order;shift];
        bank(:,k) = grid(mix(1:totaldim));
    end
end
end