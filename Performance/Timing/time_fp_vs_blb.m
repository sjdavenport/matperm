nsubj_vec = 10.^(2:6);
nvox = 1;
spfn = @(n) randn(n, nvox)';
fp_timings = zeros(1, length(nsubj_vec));
blb_timings = zeros(1, length(nsubj_vec));

alpha = 0.05;

for I = 1:length(nsubj_vec)
    I
    nsubj = nsubj_vec(I);
    data = spfn(nsubj);
    tic; fastperm_mean( data, 20, 0.05, 1000, 1, 0, 0, 1); 
    fp_timings(I) = toc;
    tic;  blb( data, floor(nsubj^(0.7)), 5, 30 );
    blb_timings7(I) = toc; 
    tic;  blb( data, floor(nsubj^(0.6)), 5, 30 );
    blb_timings6(I) = toc; 
end

%% 
save('./fp_vs_blb_timings', 'fp_timings', 'blb_timings6', 'blb_timings7', 'nsubj_vec' )

%%
clf
mploc = 'C:\Users\12SDa\davenpor\davenpor\Toolboxes\matperm\';

load([mploc, 'Performance\fp_vs_blb_timings'])
plot(log10(fp_timings))
hold on
plot(log10(blb_timings6))
plot(log10(blb_timings7))
xticklabels(nsubj_vec);
ylabel('log_{10}(Seconds)')
xlabel('Number of Subjects')
legend('Cheetah', 'BLB - gamma = 0.6', 'BLB - gamma = 0.7', 'Location', 'NorthWest')
title('Comparing the log timing') 

saveloc = [mploc, '\Performance\Figures\'];
export_fig([saveloc, 'time_fp_blb_log.pdf'])

%%
clf
plot(fp_timings)
hold on
plot(blb_timings6)
plot(blb_timings7)
xticklabels(nsubj_vec);
ylabel('Seconds')
xlabel('Number of Subjects')
legend('Cheetah', 'BLB - gamma = 0.6', 'BLB - gamma = 0.7', 'Location', 'NorthWest')
title('Comparing the total timing') 
% Note that this includes calculation of the test-statistics which is of
% course an O(n) operation. would be good to include an illustration of the
% fact that once the blocks have been calculated 

saveloc = [mploc, '\Performance\Figures\'];
export_fig([saveloc, 'time_fp_blb_usual.pdf'])

%% Time vs power
blb_timings = zeros(1, length(nsubj_vec));
nvox = 1;
spfn = @(n) randn(n, nvox)' + 0.05;
nsubj = 1000;
% Calculate timings
subset_size = 5:5:25;
store_power = zeros(1, length(subset_size));
store_times = zeros(1, length(subset_size));
for I = 1:length(subset_size)
    tic;  blb( data, floor(nsubj^(0.7)), 5, 30 ); 
    store_times(I) = toc;
    store_power(I) = rc_blb(  spfn, nsubj, 100, 5, 30, niters, alpha );
end

%%
n = 10000; nvox = 1000; data = randn(nvox,n);
subset_size = floor(n^(0.6)); n_subsamples = 50; n_mc = 100;

tic; a = blb( data, subset_size, n_subsamples, n_mc ); 
blb_time = toc;
tic; [~,c] = fastperm_mean( data, 15, 0.05, n_subsamples*n_mc ); 
fp_time = toc;

fprintf('BLB: %2.2f, FP:  %2.2f\n', blb_time, fp_time)

%%
tic; a = blb( data, subset_size, 5, 30 ); 
blb_time = toc;

%%
n = 10000; data = randn(10000,n);
subset_size = floor(n^(0.6)); n_subsamples = 50; n_mc = 100;
tic; a = blb( data, subset_size, n_subsamples, n_mc ); toc
tic; [~,c] = fastperm_mean( data, 15, 0.05, 150 ); toc
