%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the blocklmtstat function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsubj = 10; p = 4; nvox = 49;
data = randn(nvox,nsubj);
design_matrix = randn(nsubj, p);
nblocks = 11;
[ block_xY, block_sos ] = block_lm_summary_stats( data, design_matrix, nblocks );
[ betahat, sigmahat ] = blocklmtstat( block_xY, block_sos, design_matrix);
[ ~, ~, ~, betahat2, sigmahat2 ] = ...
      contrast_tstats_noerrorchecking(data, design_matrix, ones(1,p));
sum(betahat(:) - betahat2(:))
sum(sigmahat(:) - sigmahat2(:))