function [ tstat_array, Cbetahat ] = bs2tstat( betahat, sigmahat, ... 
                                                contrast_matrix, xtx_inv )
% bs2tstat( betahat, sigmahat, contrast_matrix, xtx_inv ) calculates
% voxelwise t-statistics from summary statistics of the data
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  betahat
%  sigmahat
%  contrast_matrix
%  xtx_inv 
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% nvox = 100; nsubj = 30; categ = zeros(1, nsubj);
% X = group_design(categ); C = 1; data = randn(100, nsubj); nblocks = 10;
% [ block_xY, block_sos ] = block_lm_summary_stats( data, X, nblocks );
% [ betahat, sigmahat ] = blocklmtstat( block_xY, block_sos, X);
% xtx_inv = inv(X'*X);
% tstat = bs2tstat( betahat, sigmahat, C, xtx_inv );
% tstat_ctn = contrast_tstats_noerrorchecking(data, X, C);
% sum(tstat(:) - tstat_ctn(:))
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
n_contrasts = size(contrast_matrix,1);

%%  Main Function Loop
%--------------------------------------------------------------------------
% Cbetahat = (contrast_matrix*betahat')';
Cbetahat = betahat*contrast_matrix';
tstat_array = Cbetahat./sigmahat;

% Scale by the scaling constants to ensure variance 1
for l = 1:n_contrasts
    scaling_constant = (contrast_matrix(l,:) * xtx_inv * contrast_matrix(l,:)').^(1/2);
	tstat_array(:,l) = tstat_array(:,l)/scaling_constant;
end

end

