%% Calculating the FPR for the fast permutation methods
%% Single voxel - one sample - fpr
spfn = @(n) randn(n, 1)'; nsubj = 1000;
do_mean = 1;
fpr = rc_fastperm( spfn, nsubj, 25, 1000, do_mean, 0)
fpr = rc_fastperm( spfn, nsubj, 25, 1000, do_mean, 1)


%%
fpr = rc_fastperm( spfn, nsubj, 25, 1000, do_mean, 0, 1)

%%
fpr_blb = rc_blb(  spfn, nsubj, floor(nsubj^(0.7)), 5, 30 )

%% BLB with low b doesn't control the FPR (and is already slower than the cheetah approach)
spfn = @(n) randn(n, 1)'; nsubj = 1000;
fpr_blb = rc_blb(  spfn, nsubj, floor(nsubj^(0.5)), 5, 30 )
rc_blb(  spfn, 10000, 25, 5, 30 )

%% Compare BLB and fast perm
spfn = @(n) randn(n, 1)'; nsubj = 1000;
n_subsets = 1;
fpr_blb = rc_blb(  spfn, nsubj, 25, n_subsets, 30 )  % Inflated error rate
rc_fastperm( spfn, nsubj, 25, 30, 1, 1)     % Controls the error rate

%% Single voxel - one sample - power
spfn = @(n) randn(n, 1)' + 0.01; nsubj = 100000;
do_mean = 1;
fpr_fp = rc_fastperm( spfn, nsubj, 25, 1000, do_mean, 0)
fpr_fp = rc_fastperm( spfn, nsubj, 100, 1000, do_mean, 0)
fpr_blb = rc_blb(  spfn, nsubj, floor(nsubj^(0.5)), 5, 30 )
fpr_blb = rc_blb(  spfn, nsubj, 100, 1, 30 )
%% Single voxel -orig
fpr_orig = rc_fastperm( spfn, nsubj, 1000, 1000)

%% Single voxel - one sample - power
spfn = @(n) randn(n, 1)' + 0.05; nsubj = 1000;
fpr_fp = rc_fastperm( spfn, nsubj, 25, 1000, 0, 1, 10000);
fpr_orig = rc_fastperm( spfn, nsubj, 1000, 1000, 0, 0, 10000);

%% Single voxel - one sample - fpr
nvox = 2;
spfn = @(n) randn(n, nvox)'; nsubj = 1000;
do_mean = 1;
fpr = rc_fastperm( spfn, nsubj, 25, 1000, do_mean)
fpr_blb = rc_blb(  spfn, nsubj, floor(nsubj^(0.7)), 5, 30 )

%% Single voxel - one sample - power
nvox = 2;
spfn = @(n) randn(n, 100)' + 0.05; nsubj = 1000;
do_mean = 1;
fpr_fp = rc_fastperm( spfn, nsubj, 25, 1000, do_mean)
fpr_blb = rc_blb(  spfn, nsubj, floor(nsubj^(0.7)), 5, 30 )
