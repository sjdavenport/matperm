%%
spfn = @(n) randn(n, 1)' + 0.05;
nsubj = 100; nblocks = nsubj;
niters = 1000;
nperm_set = [100, 200, 500, 1000, 2000, 3000, 4000, 5000];
perm_power = zeros(1, length(nperm_set));
for I = 1:length(nperm_set)
    perm_power(I) = rc_fastperm( spfn, nsubj, nblocks, 1000, 0, 0, niters);
end

%% 
plot(perm_power);
xticks(1:length(nperm_set));
xticklabels(nperm_set);
%%
std_error_intervals = bernstd( fpr_fp(1,:), niters, 0.95 );
hold on
h(2) = plot(std_error_intervals(1,:)', '--', 'color', 'blue');
h(3) = plot(std_error_intervals(2,:)', '--', 'color', 'blue')