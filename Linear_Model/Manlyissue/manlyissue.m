nsubj = 100;

e = randn(nsubj, 1);
beta1 = 100;
beta2 = 0;
beta = [beta1;beta2];

contrast_vector = [0,1];
n_params = length(c);

Xcol1 = [ones(nsubj/2, 1); zeros(nsubj/2, 1)];
Xcol2 = [zeros(nsubj/2, 1); ones(nsubj/2, 1)];
X = [Xcol1, Xcol2];
Y = X*beta + e;

xtx_inv = inv(X'*X);
betahat = xtx_inv * X' * Y;
rfmate = eye(nsubj) - X * xtx_inv * X';
residuals = rfmate * Y;
sigmahat = (sum(residuals.^2)*(1/(nsubj-n_params))).^(1/2);

cbetahat = betahat(1) - betahat(2);
scaling_constant = (contrast_vector * xtx_inv * contrast_vector').^(1/2);
tstatistic = cbetahat/sigmahat/scaling_constant;

%% Manly permutation step
nperm = 
for I = 