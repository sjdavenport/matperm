function [ threshold, vec_of_maxima, permuted_tstat_store ] = ...
        fastperm( data, nblocks, alpha, nperm, store_perms, init_randomize)
% FASTPERM( data, nblocks, alpha, nperm, store_perms, init_randomize)
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%    data      a Dim by nsubj array of the data
%    nblocks   a positive integer giving the number of blocks with which to
%              partition the data into
% Optional
%  alpha       a number in (0,1) giving the error rate to control at.
%              Default is 0.05.
%  nperm       a positive integer specifying the number of permutations to
%              use. Default is nperm = 5000;
%  store_perms   0/1, specifying whether to store the permuted images.
%                Default is 0, i.e. not to.
%  init_randomize  0/1 whether to perform an initial randomization of the
%                  data before partioning into blocks. Default recommended
%                  is to do so i.e. setting it to be 1.
%--------------------------------------------------------------------------
% OUTPUT
%  threshold
%  vec_of_maxima
%  permuted_tstat_store
%--------------------------------------------------------------------------
% EXAMPLES
% % Signle voxel example
% nsubj = 100; nvox = 50; nblocks = 10; alpha = 0.05; nperm = 1000;
% data = normrnd(0,1,nvox,nsubj);
% [~, threshold_orig] = perm_thresh(data, 'T', alpha, NaN, NaN, 0, nperm)
% threshold_fast_perm = fast_perm( data, nblocks, alpha, nperm )
%
% % Larger Image example
% dim = [10,10]; nsubj = 1000;
% data = noisegen(dim, nsubj, 2, 0);
% tic; [~, threshold_orig_perm] = perm_thresh(data, 'T', alpha, NaN, NaN, 0, nperm); toc
% tic; threshold_fast_perm_10blocks = fast_perm( data, 10, alpha, nperm ); toc
% tic; threshold_fast_perm_30blocks = fast_perm( data, 30, alpha, nperm ); toc
% tic; threshold_fast_perm_50blocks = fast_perm( data, 50, alpha, nperm ); toc
% threshold_orig_perm
% threshold_fast_perm_10blocks
% threshold_fast_perm_30blocks
% threshold_fast_perm_50blocks
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check optional input
%--------------------------------------------------------------------------
if ~exist('nperm', 'var')
    nperm = 5000;
end

if ~exist('alpha', 'var')
    alpha = 0.05;
end

if ~exist('store_perms', 'var')
    store_perms = 0;
end

if ~exist('init_randomize', 'var')
    init_randomize = 1;
end

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
s_data = size(data);
nsubj = s_data(end);
dim = s_data(1:end-1);
D = length(dim);

% Initialize the vector of permuted t-statistics
if store_perms == 1
    permuted_tstat_store = zeros([dim, nperm]);
else
    % Set to NaN if not storing it
    permuted_tstat_store = NaN;
end

% Initialize the max_vector
vec_of_maxima = zeros(1, nperm);

% Create a variable index
variable_index = repmat( {':'}, 1, D );

%%  Main Function Loop
%--------------------------------------------------------------------------
% Calculate the sum and the mean squared within each of the blocks
[ block_sum, block_sos ] = block_summary_stats( data, nblocks );

% Compute the original t-statistic and set that to be the first permutation
original_tstat = block_tstat(block_sum, block_sos, nsubj);

% Compute the maximal t-statistic
vec_of_maxima(1) = max(original_tstat(:));

% Initial randomization step (done to ensure that the blocks don't contain
% consistent signal
if init_randomize
    % Generate random +- 1s
    random_bern_init = 2*(binornd(1,0.5, nsubj, 1)-1/2);
    random_sample_negative_init = find(random_bern_init < 0);
    
    % Multiply the data by +- 1
    data(variable_index{:}, random_sample_negative_init) = ...
        -data(variable_index{:}, random_sample_negative_init);
    
    % Calculate the sum and the mean squared within each of the blocks of
    % the randomized data. Note this only needs to be done if you perform 
    % an initial randomization because otherwise data is unchanged so 
    % block_sum and block_sos do not change
    [ block_sum, block_sos ] = block_summary_stats( data, nblocks );
end

% Store the first entry as the original t-statistic
if store_perms == 1
    permuted_tstat_store(variable_index{:}, 1) = original_tstat;
end

% Compute bernoulli random variables for the sign flipping
random_berns = 2*(binornd(1,0.5, nblocks, nperm )-1/2);

for I = 2:nperm
    %     modul(iter, 100);
    random_berns_for_iter = random_berns(:, I);
    random_sample_negative = find(random_berns_for_iter < 0);
    
    % Note for the python implementation need to use .copy() here!
    permuted_block_sum = block_sum;
    % Note there is no need to permute the block sum of squares as
    % minus signs are cancelled out by the square!
    
    % Calculate the permuted block sum
    permuted_block_sum(variable_index{:}, random_sample_negative) = ...
        -permuted_block_sum(variable_index{:}, random_sample_negative);
    
    % Calculate the permuted t-statistics
    permuted_tstat = block_tstat(permuted_block_sum, block_sos, nsubj);
    
    % Compute the maximal permuted t-statistic
    vec_of_maxima(I) = max(permuted_tstat(:));
    
    % Store the permuted t-statistic image (if desired)
    if store_perms == 1
        permuted_tstat_store(variable_index{:}, I) = permuted_tstat;
    end
end

% Calculate the alpha level FWER threshold
threshold = prctile(vec_of_maxima, 100*(1-alpha) );

end

