spfn = @(n) randn(n, 1)' + 0.1; nsubj = 1000;
remove_data_mean = 0;
rc_blb(  spfn, nsubj, floor(nsubj^(0.7)), subset_size(I), 30, niters, alpha, remove_data_mean );

rc_blb(  spfn, nsubj, floor(nsubj^(0.7)), subset_size(I), 30, niters, alpha, remove_data_mean );
