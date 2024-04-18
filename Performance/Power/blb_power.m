%% Calculating the power for the fast permutation methods
subset_size = [5,10,20,50,100,250,1000];
nsubj = 1000;
niters = 1000;
nvox = 1;

fpr_blb = zeros(1, length(subset_size));

nsubsamples = 5;
nmc = 30;

for I = 1:length(subset_size)
    I
    spfn = @(n) randn(n, nvox)' + 0.05;
    fpr_blb(I) = rc_blb(  spfn, nsubj, subset_size(I), nsubsamples, nmc, niters );
end

mploc = 'C:\Users\12SDa\davenpor\davenpor\Toolboxes\matperm\';

save([mploc,'Performance/Poer/blb_power_vs_subsetsize_mean_nsubj_', num2str(nsubj), ...
            '_nsubsamples_', num2str(nsubsamples), '_nmc_', num2str(nmc), ...,
                '_nvox_', num2str(nvox)], 'fpr_blb', 'niters', 'subset_size')
            
%%
rc_fastperm( spfn, nsubj, 100, 1000, 1, 0, 0, niters)

%% Calculating the power for the fast permutation methods
subset_size = [5,10,20,50,100,250,1000];
nsubj = 1000;
nperm = 1000;
niters = 1000;
nvox = 100;

fpr_blb = zeros(1, length(subset_size));
nsubsamples = 20;
nmc = 30;

dist_type = 'L';

effect_size = 0;

for I = 1:length(subset_size)
    if dist_type == 'L'
        spfn = @(n) wfield([n,1], 1, 'L', 3).field' + effect_size;
    elseif dist_type == 'G'
        spfn = @(n) randn(n, nvox)' + effect_size;
    end
    fpr_blb(I) = rc_blb(  spfn, nsubj, subset_size(I), nsubsamples, nmc, niters );
end

mploc = 'C:\Users\12SDa\davenpor\davenpor\Toolboxes\matperm\';

save([mploc,'Performance/FPR/blb_fpr_vs_subsetsize_mean_dist_', ...
                                dist_type, '_nsubj_', num2str(nsubj), ...
            '_nsubsamples_', num2str(nsubsamples), '_nmc_', num2str(nmc), ...,
                '_nvox_', num2str(nvox)], 'fpr_blb', 'niters', 'subset_size')

%% Plotting FPR
mploc = 'C:\Users\12SDa\davenpor\davenpor\Toolboxes\matperm\';

nsubj = 1000;
nsubsamples = 20;
nmc = 30;
nvox = 100;
dist_type = 'L';

load([mploc, 'Performance\FPR\blb_fpr_vs_subsetsize_mean_dist_', ...
           dist_type, '_nsubj_', num2str(nsubj), '_nsubsamples_', num2str(nsubsamples), ...
                    '_nmc_', num2str(nmc), '_nvox_', num2str(nvox)])   
h(1) = plot(fpr_blb(1,:), 'color', 'blue');
std_error_intervals = bernstd( 0.05*ones(1, size(fpr_blb,2)), niters, 0.95 );
hold on
h(2) = plot(std_error_intervals(1,:)', '--', 'color', 'black');
h(3) = plot(std_error_intervals(2,:)', '--', 'color', 'black');
ylabel('FPR')
title(['BLB FPR vs subsetsize for s = ', num2str(nsubsamples), ', nmc = ', num2str(nmc)])
xticklabels(subset_size)
xlabel('Subset size')

%%
rc_blb( spfn, nsubj, subset_size(I), nsubsamples, nmc, niters );
            