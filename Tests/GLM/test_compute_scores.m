%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the compute_scores function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
nvoxels = 1000;
nsubj = 500;
nparameters = 5;  % dimension of each observation
Z = ones(nsubj, 1);
X = randn(nsubj, nparameters-1);
ZX_design = [Z,X];
gamma = randn(nvoxels, nparameters);
p = sigmoid(ZX_design * gamma');
y = binornd(1, p)';
scores = compute_scores_mv(y, Z, X, 'Binomial', 'logit');
scores_firstcol = compute_scores_mv(y(1, :), Z, X, 'Binomial', 'logit');

%% %% 3D Examples
%% Simple 3D example