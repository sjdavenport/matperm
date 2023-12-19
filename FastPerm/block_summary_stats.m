function [ block_sum, block_sos ] = block_summary_stats( data, nblocks )
% BLOCK_SUMMARY_STATS( data, nblocks ) computes the one-sample
% block summary statistics of the data for use in the implementation of
% fast permutation techniques.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  data     a matrix of size [dim, nsubj]
%  nblocks  the number of blocks to divide the data into
%--------------------------------------------------------------------------
% OUTPUT
%  block_sum   an array of size [dim, nblocks] such that block_sum(..., I) 
%              is the voxelwise sum within the Ith block
%  block_sos   an array of size [dim, nblocks] such that block_sos(..., I) 
%       is the voxelwise sum of the squares of the entries of the Ith block
%--------------------------------------------------------------------------
% EXAMPLES
% dim = [2,2]; nsubj = 100;
% randn([dim, nsubj])
% tstat_orig = mvtstat(data, dim)
% [block_sum, block_sos] = block_summary_stats( data, 10 );
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
s_data = size(data);
nsubj = s_data(end);
dim = s_data(1:end-1);
D = length(dim);

% Compute the number of subjects per block
nsubj_per_block = floor(nsubj/nblocks);
if nsubj_per_block < 1
    warning('There is less than one subject per block so defaulting to a block size of 1 subject')
    nsubj_per_block = 1;
    nblocks = nsubj;
end

% Initialize the arrays to store the block sum and sum of squares
block_sum = zeros([dim, nblocks]);
block_sos = zeros([dim, nblocks]); % sos = sum of squares

%%  Main Function Loop
%--------------------------------------------------------------------------
% Create a variable index
variable_index = repmat( {':'}, 1, D );

% Initialize the indices for the blocks
block_subject_indices = (1:nsubj_per_block) - nsubj_per_block;
for I = 1:nblocks
    block_subject_indices = block_subject_indices + nsubj_per_block;
    
    % The last block may have less subjects if the blocks are not split
    % evenly
    if I == nblocks
        block_subject_indices = block_subject_indices(block_subject_indices <= nsubj);
    end
    
    % Assign the block means
    block_sum(variable_index{:}, I) = sum(data(variable_index{:}, block_subject_indices),D+1);
    
    % Assign the block sum of squares
    block_sos(variable_index{:}, I) = sum(data(variable_index{:}, block_subject_indices).^2, D+1);
end

end

