spfn = @(n) randn(n, 1)' + 0.1; nsubj = 1000;

subset_size = 5:5:20;
alpha = 0.05;
niters = 1000;
blb_power_store = zeros(1,length(subset_size));
for I = 1:length(subset_size)
    I
    blb_power_store(I) = rc_blb(  spfn, nsubj, floor(nsubj^(0.7)), subset_size(I), 30, niters, alpha );
end