% Create y and X matrices
y = [1; 0; 1; 0; 1]';
X = [0.74, -0.042; -1.48, -0.12; 2.17, -0.55; -0.20, -1.06; -1.77, -0.30];

% Create Z matrix
Z = ones(5, 1);

% Create XZ_design by concatenating Z and X
XZ_design = [Z, X];

% Define the family (assuming 'binomial' corresponds to a binomial distribution)
distbn = 'Binomial';
linkfn = 'logit';

% Fit the Generalized Linear Model (GLM)
[betahat, fitted_values,  perfect_separations, pvalues, std_errors] = ...
                            glm_seq(y, XZ_design, distbn, linkfn);

disp(betahat')
disp(fitted_values')
disp(pvalues')
disp(std_errors')

%% Dhat and Vhat
% Calculate Dhat and Vhat
[Dhat, Vhat] = get_par_expo_fam(fitted_values', 'binomial', 'logit');

% Display Dhat and Vhat
disp(Dhat);
disp(Vhat);

%% Computing scores
% Fit the Generalized Linear Model (GLM)
[betahat0, fitted_values0] = glm_seq(y, Z, distbn, linkfn);

disp(betahat0)
disp(fitted_values0)

scores = compute_scores(y, Z, X, fitted_values0', 'Binomial', 'logit')

%%
scores_mv = compute_scores_mv(y, Z, X, 'Binomial', 'logit');
disp(squeeze(scores_mv))