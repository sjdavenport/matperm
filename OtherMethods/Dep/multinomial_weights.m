%% Full generation
subset_size = 500;
nsubj = 10000;
weights = mnrnd(nsubj, (1/subset_size)*ones(1,subset_size));

%% Subsampling (crappy)
subset_size = 500;
nsubj = 10000;
weights = (nsubj/subset_size)*mnrnd(subset_size, (1/subset_size)*ones(1,subset_size));

%% Poisson

