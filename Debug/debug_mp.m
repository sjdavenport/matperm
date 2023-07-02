nsubj = 10; p = 4; nvox = 49;
data = randn(nvox,nsubj);
design_matrix = randn(nsubj, p);
nblocks = 10;
[ block_xY, block_sos, block_Y ] = block_lm_summary_stats( data, design_matrix, nblocks, 1 );
[ betahat, sigmahat ] = blocklmtstat_orig( block_xY, block_sos, design_matrix);
[ ~, ~, ~, betahat2, sigmahat2 ] = ...
      contrast_tstats_noerrorchecking(data, design_matrix, ones(1,p));
sum(betahat(:) - betahat2(:))
sum(sigmahat(:) - sigmahat2(:)) %NOT THE SAME NEED TO INVESTIGATE!!!

%%
% Basically blocklmtstat is designed for the residuals blocklmtstat_orig for the original data!

%%
sum(block_sos,2)
sum(data.^2,2)

%%
sum(block_Y,2)
sum(data,2)

%%
