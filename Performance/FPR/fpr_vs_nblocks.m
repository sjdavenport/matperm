%% Checking the false positive rate
% nsubj = 1000;
% numberofblocks = [5,10,20,50,100,250,1000];
% nsubj = 200;
% numberofblocks = [5,10,20,50,100,200];

nsubj = 10000;
numberofblocks = [5,10,20,50,100,250,1000];

nperm = 1000;
niters = 1000;

fpr_fp = zeros(2, length(numberofblocks));

for I = 1:length(numberofblocks)
    I
    spfn = @(n) randn(n, 1)';
    fpr_fp(1,I) = rc_fastperm( spfn, nsubj, numberofblocks(I), nperm, 1, 0, niters);
    fpr_fp(2,I) = rc_fastperm( spfn, nsubj, numberofblocks(I), nperm, 1, 1, niters);
end

save(['./fpr_vs_nblocks_', num2str(nsubj)], 'fpr_fp', 'numberofblocks')

%% Plotting FPR
mploc = 'C:\Users\12SDa\davenpor\davenpor\Toolboxes\matperm\';

numberofblocks = [5,10,20,50,100,250,1000];
niters = 1000;
nsubj = 1000;
load([mploc, 'Performance\fpr_vs_nblocks_', num2str(nsubj), '.mat'])
h(1) = plot(fpr_fp(1,:), 'color', 'blue');
std_error_intervals = bernstd( 0.05*ones(1, size(fpr_fp,2)), niters, 0.95 );
hold on
h(2) = plot(std_error_intervals(1,:)', '--', 'color', 'black');
h(3) = plot(std_error_intervals(2,:)', '--', 'color', 'black');
h(4) = plot(fpr_fp(2,:), 'color', 'red');
legend([h(1), h(4)], {'No initial randomization', 'Initial randomization'}, 'Location', 'SouthEast')
ylabel('FPR')
title('FPR vs Number of Blocks')
xticklabels(numberofblocks)
xlabel('Number of blocks')

%% Checking how things vary with subject size for a fixed number of blocks
nsubj_vec = [25,50,100, 500, 1000, 5000, 10000];
nblocks = 25;

nperm = 1000;
niters = 10000;

fpr_fp = zeros(2, length(nsubj_vec));

for I = 1:length(nsubj_vec)
    I
    spfn = @(n) randn(n, 1)';
    fpr_fp(1,I) = rc_fastperm( spfn, nsubj_vec(I), nblocks, nperm, 1, 0, niters);
    fpr_fp(2,I) = rc_fastperm( spfn, nsubj_vec(I), nblocks, nperm, 1, 1, niters);
end

save(['./fpr_vs_nsubj_nblocks', num2str(nblocks)], 'fpr_fp', 'nblocks', 'nsubj_vec', 'niters')

%% Plotting FPR
mploc = 'C:\Users\12SDa\davenpor\davenpor\Toolboxes\matperm\';

nblocks = 25;
load([mploc, 'Performance\fpr_vs_nsubj_nblocks', num2str(nblocks), '.mat'])
h(1) = plot(fpr_fp(1,:), 'color', 'blue');
std_error_intervals = bernstd( 0.05*ones(1, size(fpr_fp,2)), niters, 0.95 );
hold on
h(2) = plot(std_error_intervals(1,:)', '--', 'color', 'black');
h(3) = plot(std_error_intervals(2,:)', '--', 'color', 'black');
h(4) = plot(fpr_fp(2,:), 'color', 'red');
legend([h(1), h(4)], {'No initial randomization', 'Initial randomization'}, 'Location', 'SouthEast')
ylabel('FPR')
title('FPR vs Number of Blocks')
xticklabels(nsubj_vec)
xlabel('Number of blocks')

%% For nblocks < 15 
% If the number of permutations is similar to the total number (i.e. 2^nblocks)
% then a random sample could be a biased sample potentially??

%% Studying small numbers of blocks!
numberofblocks = 5:10;
nsubj = 1000;
nperm = 1000;
niters = 1000;

fpr_fp = zeros(2, length(numberofblocks));

for I = 1:length(numberofblocks)
    I
    spfn = @(n) randn(n, 1)';
    fpr_fp(1,I) = rc_fastperm( spfn, nsubj, numberofblocks(I), nperm, 0, 0, niters);
    fpr_fp(2,I) = rc_fastperm( spfn, nsubj, numberofblocks(I), nperm, 0, 1, niters);
end

save('./fpr_vs_smallnblocks', 'fpr_fp')

%% Plotting FPR
mploc = 'C:\Users\12SDa\davenpor\davenpor\Toolboxes\matperm\';

numberofblocks = 5:10;
niters = 1000;
nsubj = 1000;
load([mploc, 'Performance\fpr_vs_smallnblocks.mat'])
h(1) = plot(fpr_fp(1,:), 'color', 'blue');
std_error_intervals = bernstd( 0.05*ones(1, size(fpr_fp,2)), niters, 0.95 );
hold on
h(2) = plot(std_error_intervals(1,:)', '--', 'color', 'black');
h(3) = plot(std_error_intervals(2,:)', '--', 'color', 'black');
h(4) = plot(fpr_fp(2,:), 'color', 'red');
legend([h(1), h(4)], {'No initial randomization', 'Initial randomization'}, 'Location', 'SouthEast')
ylabel('FPR')
title('FPR vs Number of Blocks')
xticklabels(numberofblocks)
xlabel('Number of blocks')

%%
nsubj = 10000;
nvox = 100;
spfn = @(n) randn(n, nvox)';
rc_fastperm( spfn, nsubj, 25, 1000, 0, 1, 1000)