function [tstat, xbar, std_dev, cohensd] = block_tstat( block_sum, block_sos, nsubj )
% BLOCK_TSTAT( block_means, block_sos ) takes in the means and sum of
% squares of the different blocks/partitions and computes a one-sample
% t-statistic.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  block_means      an array which is dim by nblocks for some vector dim 
%                   such that each block_means(:, ..., :, I) is the mean of
%                   the ith block
%                   This can be calculated using block_summary_stats
%  block_sos        an array which is dim by nblocks for some vector dim 
%                   such that each block_sos(:, ..., :, I) is the sum of 
%                   squares of the ith block
%                   This can be calculated using block_lm_summary_stats
%  nsubj            the number of original subjects
%--------------------------------------------------------------------------
% OUTPUT
% tstat         the one sample t-statistic at each voxel
% xbar          the mean at each voxel
% std_dev       the standard deviation at each voxel (calculated using the
%               unbiased estimate of the variance)
%--------------------------------------------------------------------------
% EXAMPLES
% % An example to ensure the output is correct
% dim = [2,2]; nsubj = 100;
% randn([dim, nsubj])
% tstat_orig = mvtstat(data, dim)
% [block_sum, block_sos] = block_summary_stats( data, 10 );
% tstat_block = block_tstat( block_sum, block_sos, nsubj )
%
% % Compare the time for large images (note that in practice the block_summary_stats
% % will only need to be calculated once wen doing permutations 
% so we exclude it from the comparison) (takes a while to generate so many images!)
% dim = [1000,1000]; nsubj = 100;
% data = noisegen(dim, nsubj, 2, 0);
% tic
% tstat_orig = mvtstat(data, dim);
% toc
% [block_sum, block_sos] = block_summary_stats( data, 10 );
% tic
% tstat_block = block_tstat( block_sum, block_sos, nsubj );
% toc
%
% % Compare the time for large images (note that in practice the block_summary_stats
% % will only need to be calculated once wen doing permutations 
% so we exclude it from the comparison) (takes a while to generate so many images!)
% dim = [10,10]; nsubj = 1000;
% data = noisegen(dim, nsubj, 2, 0);
% tic
% tstat_orig = mvtstat(data, dim);
% toc
% [block_sum, block_sos] = block_summary_stats( data, 10 );
% tic
% tstat_block = block_tstat( block_sum, block_sos, nsubj );
% toc
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Calculate the size of the input
s_block_sum = size(block_sum);

% Calculate the dimension
D = length(s_block_sum) - 1;

%%  Main Function Loop
%--------------------------------------------------------------------------
% Compute the mean
xbar = (1/nsubj)*sum(block_sum, D + 1);

% Compute the (biased) estimate of the variance
biased_variance = (1/nsubj)*sum(block_sos, D + 1) - xbar.^2;

% Compute the standard deviation using the unbiased estimator of the
% variance
std_dev = sqrt((nsubj/(nsubj-1))*biased_variance);

% Compute Cohen's d
cohensd = xbar./std_dev;

% Compute the t-statistic
tstat = cohensd*sqrt(nsubj);

end

