numberofblocks = [5,10,20,50,100,250,1000];
nsubj = 1000;
nperm = 1000;
niters = 1000;
nvox = 100;

fpr_fp = cell(2,2);

domeanperm = 1;

% Initialize the storage vectors
for init_randomize = 0:1
    for demean = 0:1
        fpr_fp{init_randomize+1, demean+1} = zeros(1, length(numberofblocks));
    end
end

for I = 1:length(numberofblocks)
    I
    spfn = @(n) randn(n, nvox)' + 0.05;
    store_index = 0;
    for init_randomize = 0:1
        for demean = 0:1
            fpr_fp{init_randomize+1, demean+1}(I) = rc_fastperm( spfn, nsubj, numberofblocks(I), nperm, domeanperm, init_randomize, demean, niters);
        end
    end
end

mploc = 'C:\Users\12SDa\davenpor\davenpor\Toolboxes\matperm\';
if domeanperm == 1
    save([mploc,'Performance/Power/power_vs_nblocks_mean_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox)], 'fpr_fp', 'niters', 'numberofblocks')
else
    save([mploc,'Performance/Power/power_vs_nblocks_tstat_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox)], 'fpr_fp', 'niters', 'numberofblocks')
end