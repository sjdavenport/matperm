%% Calculating the variance in the decided thresholds for the fast permutation methods
numberofblocks = [5,10,20,50,100,250,1000];
nsubj = 1000;
nperm = 1000;
niters = 100;

threshold_var = zeros(1, length(numberofblocks));
alpha = 0.05;
spfn = @(n) randn(n, 1)';
demean = 0;
show_loader = 0;
for I = 1:length(numberofblocks)
    data = spfn(nsubj);
    thresholds = zeros(1, niters);
    for J = 1:niters
        J
        thresholds(J) = fastperm_mean( data, numberofblocks(I), alpha, nperm, demean, show_loader, 0, 1);
    end
    threshold_var(I) = var(thresholds);
end

save('./threshvar_vs_nblocks', 'threshold_var')

%%
% load('./threshvar_vs_nblocks', 'threshold_var')
plot(threshold_var)
xticks(1:length(threshold_var));
xticklabels(numberofblocks);

%% Varying the number of permutations
nsubj = 1000;
niters = 100;
nblocks = 100;
nperm_vec = [100, 500, 1000,2000, 5000, 10000, 100000];

threshold_var = zeros(1, length(nperm_vec));
alpha = 0.05;
spfn = @(n) randn(n, 1)';
demean = 0;
show_loader = 0;
for I = 1:length(nperm_vec)
    data = spfn(nsubj);
    thresholds = zeros(1, niters);
    for J = 1:niters
        J
        thresholds(J) = fastperm_mean( data, nblocks, alpha, nperm_vec(I), demean, show_loader, 0, 1);
    end
    threshold_var(I) = var(thresholds);
end

save('./threshvar_vs_nperm', 'threshold_var')

%%
plot(threshold_var)
ylim([0,0.07])
xticks(1:length(threshold_var));
xticklabels(nperm_vec);

%% Varying the number of subjects
nsubj_vec = [100,1000,5000,10000,50000,100000];
nperm = 1000;
nblocks = 100;
niters = 100;

threshold_var = zeros(1, length(nsubj_vec));
alpha = 0.05;
spfn = @(n) randn(n, 1)';
demean = 0;
show_loader = 0;
for I = 1:length(nsubj_vec)
    data = spfn(nsubj);
    thresholds = zeros(1, niters);
    for J = 1:niters
        J
        thresholds(J) = fastperm_mean( data, nblocks, alpha, nperm, demean, show_loader, 0, 1);
    end
    threshold_var(I) = var(thresholds);
end

save('./threshvar_vs_nsubj', 'threshold_var')

%%
plot(threshold_var)
xticks(1:length(threshold_var));
xticklabels(nsubj_vec);
ylim([0, max(threshold_var)*1.2])
xlabel('nsubj')


%%
nvox = 100; nsubj = 1000;
data = normrnd(0,1,nvox,nsubj);
%%
nvox = 100; nsubj = 1000; alpha = 0.05;
out = blb( data, floor(nsubj^(0.7)), 5, 30 );
blb_thresh = prctile(out, 100*(1-alpha) );
nblocks = 25; nperm = 1000; show_loader = 0; store_perms = 0; init_randomize = 0;
fp_thresh = fastperm_mean( data, nblocks, alpha, nperm, 1, show_loader, store_perms, init_randomize);

fprintf('BLB: %1.2f, FP: %1.2f\n', blb_thresh, fp_thresh)

%%
nvox = 100; nsubj = 10000;
data = normrnd(0,1,nvox,nsubj);
nblocks = 100;
%%
fp_thresh = fastperm_mean( data, 200, alpha, 10000, 1, show_loader, 0, 1)
