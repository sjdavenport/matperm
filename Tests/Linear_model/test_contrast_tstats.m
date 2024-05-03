%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the contrast_tstats function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

age = randi([40,60], nsubj, 1);
age = age-mean(age);
sex = rand(nsubj, 1) > 0.6;
sex = sex - mean(sex);
Z = [ones(nsubj,1), sex, age];
X = randn(nsubj,1);
design_matrix = [X,Z];
nvox = 1000;
X_coeff = 2;
data = 5 + 0.2*sex' + 0.12*age' + X_coeff*X' + randn(nvox, nsubj);
[tstat_array, residuals, Cbetahat, betahat, sigmahat ] = contrast_tstats( data, design_matrix, [1,0,0,0] );
mean(betahat)
histogram(tstat_array)

%%
X_coeff = 0;
data = 5 + 0.2*sex' + 0.12*age' + X_coeff*X' + randn(nvox, nsubj);
[tstat_array, residuals, Cbetahat, betahat, sigmahat ] = contrast_tstats( data, design_matrix, [1,0,0,0] );
mean(betahat)
histogram(tstat_array)