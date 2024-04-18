function [ threshold, vec_of_maxima, permuted_tstat_store ] = ... 
      fastlmperm( data, design_matrix, contrast_matrix, nblocks, alpha,...
                                        nperm, store_perms, init_randomize)
% FASTLMPERM(  data, design_matrix, contrast_matrix, nblocks, alpha, nperm, ...
% store_perms, init_randomize)
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  data      an array of size [nvox, nblocks]
%  design_matrix    an nsubj by p matrix giving the covariates (p being the 
%                   number of parameters
%  contrast_matrix: an L by p matrix corresponding to the contrast matrix, 
%                  such that which each row is a contrast vector (where L 
%                  is the number of constrasts)
%  nblocks   a positive integer giving the number of blocks with which to
%            partition the data into
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
% nsubj = 10; p = 4; nvox = 49;
% data = randn(nvox,nsubj);
% design_matrix = randn(nsubj, p);
% nblocks = 10; contrast_matrix = ones(1,p);
% [ threshold, vec_of_maxima, permuted_tstat_store ] = ...
%              fastlmperm( data, design_matrix, contrast_matrix, nblocks );
% 
% % One sample comparison to perm_thresh
% nsubj = 10; nvox = 49;
% data = randn(nvox,nsubj);
% design_matrix = ones(nsubj,1); contrast_matrix = 1;
% threshold_flp = fastlmperm( data, design_matrix, contrast_matrix, nblocks, 0.05, 1000 )
% [~, threshold_pm] = perm_thresh( data, 'T' )
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
nvox = s_data(1);
xtx_inv = inv(design_matrix'*design_matrix);

% Initialize the vector of permuted t-statistics
if store_perms == 1
    permuted_tstat_store = zeros([nvox, nperm]);
else
    % Set to NaN if not storing it
    permuted_tstat_store = NaN;
end

% Initialize the max_vector
vec_of_maxima = zeros(1, nperm);

% Performing error checking on the linear model input
contrast_tstats_errorchecking( data, design_matrix, contrast_matrix );

%%  Main Function Loop
%--------------------------------------------------------------------------
% Calculate the sum and the mean squared within each of the blocks
doblockY = 1;
[ block_xY, block_sos, block_Y ] = block_lm_summary_stats( data, ... 
                                        design_matrix, nblocks, doblockY );

% Compute betahat and sigmahat
error('I think the line below should be blocklmtstat_orig and without the block Y perhaps??, need to investigate!')
[ betahat, sigmahat ] = blocklmtstat( block_xY, block_sos, block_Y, design_matrix);

% Compute the t-statistics on the original data
original_tstat = bs2tstat( betahat, sigmahat, contrast_matrix, xtx_inv );

% Compute the maximal t-statistic
vec_of_maxima(1) = max(original_tstat(:));

% Calculate the residual forming matrix
rf_matrix = eye(nsubj) - design_matrix * xtx_inv * design_matrix';

% Compute the residuals (I-P)Y
residuals = data*rf_matrix';

% Initial randomization step (done to ensure that the blocks don't contain
% consistent signal
if init_randomize
    % Generate random +- 1s
    random_bern_init = 2*(binornd(1,0.5, nsubj, 1)-1/2);
    random_sample_negative_init = find(random_bern_init < 0);
    
    % Multiply the residuals by +- 1
    residuals(:, random_sample_negative_init) = -residuals(:, random_sample_negative_init);
    
    % Calculate the sum and the mean squared within each of the blocks
    [ block_xY, block_sos ] = block_lm_summary_stats( residuals, design_matrix, nblocks );
end

% Store the first entry as the original t-statistic
if store_perms == 1
    permuted_tstat_store(:, 1) = original_tstat;
end

% Compute bernoulli random variables for the sign flipping
random_berns = 2*(binornd(1,0.5, nblocks, nperm )-1/2);

for I = 2:nperm
    %     modul(iter, 100);
    random_berns_for_iter = random_berns(:, I);
    random_sample_negative = find(random_berns_for_iter < 0);
    
    % Note for the python implementation need to use .copy() here!
    permuted_block_xY = block_xY;
    permuted_block_Y = block_Y;
    % Note there is no need to permute the block sum of squares as
    % minus signs are cancelled out by the square!
    
    % Calculate the permuted block sum
    permuted_block_xY(:, :, random_sample_negative) = ...
                        -permuted_block_xY(:, :, random_sample_negative);
                    
    permuted_block_Y(:, random_sample_negative) = ...
                             -permuted_block_Y(:, random_sample_negative);
    
    % Calculate the permuted t-statistics
    [ betahat_perm, sigmahat_perm ] = ...
         blocklmtstat(permuted_block_xY, block_sos, permuted_block_Y, ...
                                                            design_matrix);
    permuted_tstat = bs2tstat( betahat_perm, sigmahat_perm,...
                                                contrast_matrix, xtx_inv );

    % Compute the maximal permuted t-statistic
    vec_of_maxima(I) = max(permuted_tstat(:));
    
    % Store the permuted t-statistic image (if desired)
    if store_perms == 1
        permuted_tstat_store(:, I) = permuted_tstat;
    end
end

% Calculate the alpha level FWER threshold
threshold = prctile(vec_of_maxima, 100*(1-alpha) );

end

