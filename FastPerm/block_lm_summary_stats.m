function [ block_xY, block_sos, block_Y ] = ...
                  block_lm_summary_stats( data, design, nblocks, doblockY )
% BLOCK_LM_SUMMARY_STATS( lat_data, design, contrast_matrix, nblocks )
% computes X^TY and Y.^2 in blocks for input into fast permutation
% functions
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% data     a matrix of size [nvox, nsubj]
% design:    an N by p matrix giving the covariates (p being the number of
%            parameters
% nblocks:   the number of blocks to divide the subjects into - this should
%            be less than or equal to the number of subjects
%--------------------------------------------------------------------------
% OUTPUT
% block_xY   an array with dimension given by the vector [nvox, p, nblocks]
%            sum that block_xY(..., I) = \sum_{i in Ith block} x_i Y_i
%            for I = 1:nblocks, where x_i^T is the ith row of X (1 \leq i
%            \leq nsubj).
% block_sos   an array with dimension given by the vector [nvox, nblocks] 
%            such that block_sos(..., I) is the voxelwise sum of the squares 
%            of the entries of the Ith block ie:
%            \sum_{i in Ith block} Y_i.^2
% block_Y    an array with dimension given by the vector [nvox, nblocks]
%            sum that block_Y(:, I) = \sum_{i in Ith block} Y_i
%            for I = 1:nblocks.
%--------------------------------------------------------------------------
% EXAMPLES
% nsubj = 10; p = 4; nvox = 49;
% lat_data = randn([nvox,nsubj]);
% design = randn(nsubj, p);
% nblocks = 5;
% [ block_xY, block_sos ] = block_lm_summary_stats( lat_data, design, nblocks )
% [ betahat, sigmahat ] = blocklmtstat_orig( block_xY, block_sos, design)
% tstat = 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

if ~exist('doblockY', 'var')
    doblockY = 0;
end

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
s_data = size(data);
nsubj = s_data(end);
dim = s_data(1:end-1);
D = length(dim);

% Calculate the number of parameters
nparams = size(design,2);

% Compute the number of subjects per block
nsubj_per_block = ceil(nsubj/nblocks);
if nsubj_per_block < 1
    warning('There is less than one subject per block so defaulting to a block size of 1 subject')
    nsubj_per_block = 1;
    nblocks = nsubj;
end

% Initialize the arrays to store the block sum and sum of squares
block_xY = zeros([dim, nparams, nblocks]);
block_sos = zeros([dim, nblocks]); % sos = sum of squares
if doblockY
    block_Y = zeros([dim, nblocks]);
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% Create a variable index
variable_index = repmat( {':'}, 1, D);

% Initialize the indices for the blocks
block_subject_indices = (1:nsubj_per_block) - nsubj_per_block;

% Main loop
for I = 1:nblocks
    block_subject_indices = block_subject_indices + nsubj_per_block;
    
    % The last block may have less subjects if the blocks are not split
    % evenly
    if I == nblocks
        block_subject_indices = block_subject_indices(block_subject_indices <= nsubj);
    end
    
    % Obtain the blocks of X^TY
%     block_xY(variable_index{:}, :, I) = marray(design(block_subject_indices, :)', data(variable_index{:}, block_subject_indices));
    block_xY(variable_index{:}, :, I) = data(variable_index{:}, block_subject_indices)*design(block_subject_indices, :);

    % Assign the block sum of squares
    block_sos(variable_index{:}, I) = sum(data(variable_index{:}, block_subject_indices).^2, D+1);    
end

if doblockY
    % Initialize the indices for the blocks
    block_subject_indices = (1:nsubj_per_block) - nsubj_per_block;

    for J = 1:nblocks
        block_subject_indices = block_subject_indices + nsubj_per_block;
        
        % Calculate the sum of the Ys
        block_Y(variable_index{:}, J) = sum(data(variable_index{:}, block_subject_indices), D+1);
    end
else
    block_Y = NaN;
end

end

