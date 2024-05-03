function [ betahat, sigmahat ] = blocklmtstat( block_xY, block_sos, ...,
                                                            design_matrix)
% BLOCKLMTSTAT( block_xY, block_sos, contrast_matrix, nsubj ) takes block
% X^TY and the sum of Y squared (voxelwise) and uses these to compute the
% test-statistics.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  block_xY    an array of size [nvox, nblocks] such that
%             block_xY(..., I) is \sum_{i in Ith block} x_i Y_i
%             This can be calculated using block_lm_summary_stats
%  block_sos  an array of size [nvox, nblocks] such that block_sos(..., I)
%             is the voxelwise sum of the squares of the entries of the Ith 
%             block. This can be calculated using block_lm_summary_stats.
%  design_matrix  the n by p design matrix where n is the number of subjects
%                and p is the number of parameters in the model
%--------------------------------------------------------------------------
% OUTPUT
%  betahat    a matrix of size nvox by p giving the length p betahat vector 
%             at each voxel
%  sigmahat   a vector of length nvox giving the standard error at each
%             voxel
%--------------------------------------------------------------------------
% EXAMPLES
% nsubj = 10; p = 4; nvox = 49;
% data = randn(nvox,nsubj);
% design_matrix = randn(nsubj, p);
% nblocks = 10;
% [ block_xY, block_sos ] = block_lm_summary_stats( data, design_matrix, nblocks );
% [ betahat, sigmahat ] = blocklmtstat( block_xY, block_sos, design_matrix);
% [ ~, ~, ~, betahat2, sigmahat2 ] = ...
%       contrast_tstats_noerrorchecking(data, design_matrix, ones(1,p));
% sum(betahat(:) - betahat2(:))
% sum(sigmahat(:) - sigmahat2(:))
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Calculate the size of the input
s_block_sum = size(block_xY);

% Calculate the dimension
D = length(s_block_sum) - 2;

% Calculate the number of subjects
nsubj = size(design_matrix,1);

%%  Main Function Loop
%--------------------------------------------------------------------------
% Compute betahat
xtx = design_matrix'*design_matrix;
xtx_inv = inv(xtx);
block_xY_sum = sum(block_xY, D + 2);
betahat = block_xY_sum*xtx_inv;  %Note no need to transpose xtx_inv because it is symmetric

% Compute sigma
sum_of_sos = sum(block_sos, D + 1);
middle_term = -2*sum(block_xY_sum.*betahat, 2);
end_term = betahat*xtx;
end_term = sum(end_term.*betahat, 2);
design_rank = rank(design_matrix);
sigmahat = (1/(nsubj - design_rank))*(sum_of_sos + middle_term + end_term);
sigmahat = sigmahat.^(1/2);

end

