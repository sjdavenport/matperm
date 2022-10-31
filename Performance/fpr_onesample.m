%% Calculating the FPR for the fast permutation methods
%% Single voxel - one sample - fpr
spfn = @(n) randn(n, 1)'; nsubj = 1000;
fpr = rc_fast_perm( spfn, nsubj, 25, 1000)

%% Single voxel - one sample - power
spfn = @(n) randn(n, 1)' + 0.05; nsubj = 1000;
fpr_fp = rc_fast_perm( spfn, nsubj, 25, 1000)
fpr_orig = rc_fast_perm( spfn, nsubj, 1000, 1000)

%% Single voxel - one sample - power
spfn = @(n) randn(n, 1)' + 0.05; nsubj = 1000;
fpr_fp = rc_fast_perm( spfn, nsubj, 25, 1000, 1, 10000);
fpr_orig = rc_fast_perm( spfn, nsubj, 1000, 1000, 0, 10000);