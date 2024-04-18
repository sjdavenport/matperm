%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the glm_seq function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Single voxel
nsubj = 10000;  % number of observations
nparameters = 5;  % dimension of each observation
X = [ones(nsubj, 1),randn(nsubj, nparameters-1)];
gamma = randn(nparameters, 1);
p = sigmoid(X * gamma);
y = binornd(1, p)';
[gammahat, fitted_values, ~, ~, ~] = glm_seq(y, X, 'Binomial', 'logit');
disp(gammahat);
disp(gamma.');

%% Single voxel Perfect separation
nsubj = 50;  % number of observations
nparameters = 2;  % dimension of each observation
X = [zeros(nsubj/2,1);ones(nsubj/2,1)];
X = X - mean(X);
X = [ones(nsubj, 1), X];
gamma = [-2,-5]';
p = sigmoid(X * gamma);
y = binornd(1, p)';
[gammahat, fitted_values, ~, ~, ~] = glm_seq(y, X, 'Binomial', 'logit');
disp(gammahat);
disp(gamma.');

%% 
% Multiple voxel example
nvoxels = 2;
nsubj = 100000;
nparameters = 5;  % dimension of each observation
X = [ones(nsubj, 1),randn(nsubj, nparameters-1)];
gamma = randn(nvoxels, nparameters);
p = sigmoid(X * gamma');
y = binornd(1, p)';
[gammahat, fitted_values, ~, ~, ~] = glm_seq(y, X, 'Binomial', 'logit');
disp(gammahat-gamma);

%% Many voxel example
nvoxels = 1000;
nsubj = 500;
nparameters = 5;  % dimension of each observation
X = [ones(nsubj, 1),randn(nsubj, nparameters-1)];
gamma = randn(nvoxels, nparameters);
p = sigmoid(X * gamma');
y = binornd(1, p)';
[gammahat, fitted_values, pslocs, ~, ~] = glm_seq(y, X, 'Binomial', 'logit');
plot(gamma(:), gammahat(:), '*')
hold on 
plotid