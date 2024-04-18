%% Calculating the power for the fast permutation methods
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
    save([mploc,'Performance/Power/power_vs_nblocks_mean_nsubj_', num2str(nsubj), ...
        '_nvox_', num2str(nvox),  '_nperm_', num2str(nperm)], 'fpr_fp', ...
                                        'niters', 'numberofblocks')
else
    save([mploc,'Performance/Power/power_vs_nblocks_tstat_nsubj_', num2str(nsubj), ...
        '_nvox_', num2str(nvox), '_nperm_', num2str(nperm)], 'fpr_fp', 'niters', ...
                                                        'numberofblocks')
end

%%
mploc = 'C:\Users\12SDa\davenpor\davenpor\Toolboxes\matperm\';

numberofblocks = [5,10,20,50,100,250,1000];
nsubj = 1000;
nvox = 1;
nperm = 1000;
load([mploc, 'Performance\Power\power_vs_nblocks_mean_nsubj_', ...
                        num2str(nsubj), '_nvox_', num2str(nvox), ...
                                        '_nperm_', num2str(nperm), '.mat'])
for demean = 0:1
    subplot(1,2, demean+1)
    h(1) = plot(fpr_fp{1,demean + 1}, 'color', 'blue');
    std_error_intervals = bernstd( fpr_fp{1,demean + 1}, niters, 0.95 );
    hold on
    h(2) = plot(std_error_intervals(1,:)', '--', 'color', 'blue');
    h(3) = plot(std_error_intervals(2,:)', '--', 'color', 'blue');
    h(4) = plot(fpr_fp{2,demean + 1}, 'color', 'red');
    std_error_intervals_2 = bernstd( fpr_fp{2,demean + 1}, niters, 0.95 );
    h(5) = plot(std_error_intervals_2(1,:)', '--', 'color', 'red');
    h(6) = plot(std_error_intervals_2(2,:)', '--', 'color', 'red');
    legend([h(1), h(4)], {'No initial randomization', 'Initial randomization'}, 'Location', 'SouthEast')
    ylabel('Power')
    if demean == 1
        title('Demeaned power')
    else
        title('Original power')
    end
    xticks(1:length(numberofblocks));
    xticklabels(numberofblocks)
    if nvox == 1
        ylim([0, 0.6])
    elseif nvox == 100
        ylim([0.9, 1])
    end
    xlabel('Number of blocks')
end
saveloc = [mploc, '\Performance\Figures\'];
export_fig([saveloc, 'power_vs_nblocks.pdf'], '-transparent')

%% Plotting power
mploc = 'C:\Users\12SDa\davenpor\davenpor\Toolboxes\matperm\';

numberofblocks = [5,10,20,50,100,250,1000];
niters = 10000;
load([mploc, 'Performance\Power\power_vs_nblocks_1000_orig.mat'])
h(1) = plot(fpr_fp(1,:), 'color', 'blue');
std_error_intervals = bernstd( fpr_fp(1,:), niters, 0.95 );
hold on
h(2) = plot(std_error_intervals(1,:)', '--', 'color', 'blue');
h(3) = plot(std_error_intervals(2,:)', '--', 'color', 'blue');
h(4) = plot(fpr_fp(2,:), 'color', 'red');
std_error_intervals_2 = bernstd( fpr_fp(2,:), niters, 0.95 );
h(5) = plot(std_error_intervals_2(1,:)', '--', 'color', 'red');
h(6) = plot(std_error_intervals_2(2,:)', '--', 'color', 'red');
legend([h(1), h(4)], {'No initial randomization', 'Initial randomization'}, 'Location', 'SouthEast')
ylabel('Power')
title('Power vs Number of Blocks')
xticklabels(numberofblocks)
xlabel('Number of blocks')

saveloc = [mploc, '\Performance\Figures\'];
export_fig([saveloc, 'power_vs_nblocks.pdf'], '-transparent')