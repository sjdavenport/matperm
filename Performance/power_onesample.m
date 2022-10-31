%% Calculating the power for the fast permutation methods
numberofblocks = [5,10,20,50,100,250,1000];
nsubj = 1000;
nperm = 1000;
niters = 10000;

fpr_fp = zeros(2, length(numberofblocks));

for I = 1:length(numberofblocks)
    spfn = @(n) randn(n, 1)' + 0.05;
    fpr_fp(1,I) = rc_fastperm( spfn, nsubj, numberofblocks(I), nperm, 0, niters);
    fpr_fp(2,I) = rc_fastperm( spfn, nsubj, numberofblocks(I), nperm, 1, niters);
end

save('./onesamplepower', 'fpr_fp')

%% Plotting power
numberofblocks = [5,10,20,50,100,250,1000];
niters = 10000;
load('./onesamplepower.mat')
h(1) = plot(fpr_fp(1,:), 'color', 'blue');
std_error_intervals = bernstd( fpr_fp(1,:), niters, 0.95 );
hold on
h(2) = plot(std_error_intervals(1,:)', '--', 'color', 'blue');
h(3) = plot(std_error_intervals(2,:)', '--', 'color', 'blue');
h(4) = plot(fpr_fp(2,:), 'color', 'red');
std_error_intervals_2 = bernstd( fpr_fp(2,:), niters, 0.95 );
h(5) = plot(std_error_intervals_2(1,:)', '--', 'color', 'red');
h(6) = plot(std_error_intervals_2(2,:)', '--', 'color', 'red');
legend([h(1), h(4)], {'No randomization', 'Initial randomization'}, 'Location', 'SouthEast')
ylabel('Power')
title('Power vs Number of Blocks')
xticklabels(numberofblocks)
xlabel('Number of blocks')
