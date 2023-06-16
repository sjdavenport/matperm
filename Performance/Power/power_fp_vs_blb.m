% spfn = @(n) randn(n, 1)' + 0.1; nsubj = 100;
nvox = 100;
spfn = @(n) randn(n, nvox)' + 0.01; nsubj = 10000;

do_mean = 1;
niters = 1000;
fpr_fp = rc_fastperm( spfn, nsubj, 15, 1000, do_mean, 1, niters, alpha);
fpr_blb = rc_blb(  spfn, nsubj, floor(nsubj^(0.7)), 5, 30, niters, alpha );

fprintf('BLB: %2.3f, FP:  %2.3f\n', fpr_blb, fpr_fp)

% In the nsubj = 100, nvox = 1 case BLB isn't controlling the FPR!! So this
% setting shouldn't be used to compare power!! If its not clear whether it
% will control the FPR then it seems inadvisable to use it for these sample
% sizes.

% Even for high gamma and nsubj and 

% Need to code up another measure of power for the 

%% For large image sizes BLB isn't able to estimate the  
nsubj = 50;
spfn = @(n) randn(n, 100)';
fpr_blb = rc_blb(  spfn, nsubj, floor(nsubj^(0.7)), 5, 30, niters, alpha )
