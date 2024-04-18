%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the firth_regression_seq function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Single voxel
nsubj = 10000;  % number of observations
nparameters = 5;  % dimension of each observation
X = [ones(nsubj, 1),randn(nsubj, nparameters-1)];
gamma = randn(nparameters, 1);
p = sigmoid(X * gamma);
y = binornd(1, p)';
[gammahat, fitted_values] = firth_regression_seq(y, X);
disp(gammahat');
disp(gamma');

%% Single voxel Perfect separation
nsubj = 50;  % number of observations
nparameters = 2;  % dimension of each observation
X = [zeros(nsubj/2,1);ones(nsubj/2,1)];
X = X - mean(X);
X = [ones(nsubj, 1), X];
gamma = [-2,-5]';
p = sigmoid(X * gamma);
y = binornd(1, p)';
[gammahat, fitted_values] = firth_regression_seq(y, X);
disp(gammahat');
disp(gamma.');

%% Many voxel example
nvoxels = 1000;
nsubj = 500;
nparameters = 5;  % dimension of each observation
X = [ones(nsubj, 1),randn(nsubj, nparameters-1)];
gamma = randn(nvoxels, nparameters);
p = sigmoid(X * gamma');
y = binornd(1, p)';
[gammahat, fitted_values] = firth_regression_seq(y, X);
plot(gamma(:), gammahat(:), '*')
hold on 
plotid