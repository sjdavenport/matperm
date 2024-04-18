%% Calculating the variance in the decided thresholds for the fast permutation methods
gamma = 0.1:0.1:1;
nsubj = 1000;
nperm = 1000;
niters = 100;

threshold_var = zeros(1, length(gamma));
alpha = 0.05;
spfn = @(n) randn(n, 1)';
demean = 0;
show_loader = 0;
for I = 1:length(gamma)
    I
    data = spfn(nsubj);
    thresholds = zeros(1, niters);
    subset_size = floor(nsubj^(gamma(I)));
    for J = 1:niters
        out = blb( data, subset_size, 20, 40 );
        thresholds(J) = prctile(out, 100*(1-alpha) );
    end
    threshold_var(I) = var(thresholds);
end

save('./threshvar_vs_blb_subsetsize', 'threshold_var')

%%
load('./threshvar_vs_blb_subsetsize', 'threshold_var')

plot(threshold_var)
xticks(1:length(threshold_var));
% xlim([1, length(gamma)]);
xticklabels(gamma);
