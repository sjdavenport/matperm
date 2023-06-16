function [ bootstrap_samples ] = blb( data, subset_size, n_subsamples, n_mc )
% BLB implements the bag of little bootstraps method
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  data     a matrix of size [dim, nsubj]
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% n = 10000; data = randn(1,n);
% subset_size = floor(n^(0.7)); n_subsamples = 50; n_mc = 100;
% bld_dist = blb( data, subset_size, n_subsamples, n_mc );
% boot_dist = bootstrap(data, n_subsamples*n_mc + 1);
% [~,fp_dist] = fastperm_mean( data, 12, 0.05, n_subsamples*n_mc + 1);
% subplot(3,1,1); histogram(bld_dist, 'BinWidth', 0.2); title('Bag of Little Bootstraps'); xlim([-4,4])
% subplot(3,1,2); histogram(boot_dist, 'BinWidth', 0.2);  title('Bootstraps'); xlim([-4,4])
% subplot(3,1,3); histogram(fp_dist, 'BinWidth', 0.2);  title('Fast Perm'); xlim([-4,4])
%
% n = 10000; data = randn(1000,n);
% subset_size = floor(n^(0.7)); n_subsamples = 25; n_mc = 50;
% tic; a = blb( data, subset_size, n_subsamples, n_mc ); toc
% tic; [~,c] = fastperm_mean( data, 15, 0.05, n_subsamples*n_mc ); toc
%
% tic; b = bootstrap(data, n_subsamples*n_mc); toc
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

if ~exist('remove_data_mean', 'var')
    remove_data_mean = 0;
end

%%  Check mandatory input and get important constants
%-------------------------------------------------------------------------
nsubj = size(data, 2);

%%  Main Function Loop
%--------------------------------------------------------------------------
bootstrap_samples = zeros(n_subsamples*n_mc, 1);
actual_mean = mean(data,2);

for I = 1:n_subsamples
    sample_subset = randsample(nsubj,subset_size,0);
    sample_data = data(:, sample_subset);
    if remove_data_mean == 0
        bootstrap_subsample =  sample_data - mean(sample_data, 2);
    end
    for J = 1:n_mc
        % Computing the weights below appears to be quite expensive! But
        % unavoidable!
         weights = mnrnd(nsubj, (1/subset_size)*ones(1,subset_size));
%          weights = floor((nsubj/1)*mnrnd(1, (1/subset_size)*ones(1,subset_size)));
         weighted_sample = sum(bootstrap_subsample.*weights,2)/nsubj;
         bootstrap_samples(n_mc*(I-1) + J) = max(sqrt(nsubj)*weighted_sample(:));
    end
end

bootstrap_samples = [max(sqrt(nsubj)*actual_mean(:)); bootstrap_samples];

end

