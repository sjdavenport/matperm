function scores = compute_scores(y, Z, X, fitted_values_0, family, link, score_type)
% COMPUTE_SCORES(y, Z, X, fitted_values_0, family, link, score_type)
% computes the score contributions for a generalized linear model.
%--------------------------------------------------------------------------
% INPUT
%   y (numeric): A matrix of shape (n_samples, 1) representing the dependent variable.
%   Z (numeric): An array of shape (n_samples, n_nuisance_parameters) representing the nuisance variables.
%   X (numeric): An array of shape (n_samples, n_parameters) representing the independent variables.
%   fitted_values_0 (numeric): Fitted values
%   family (char): A string specifying the distribution of the response variable. Must be one of:
%                  'Binomial', 'Gamma', 'Gaussian', 'Inverse Gaussian', 'Negative Binomial', 'Poisson'.
%   link (char): A string specifying the link function. Must be one of: 'log', 'logit', 'probit', 'cauchy', 'cloglog', 'identity', 'inverse'.
%   score_type (char, optional): A string specifying the type of scores to compute. Default is 'effective'.
%--------------------------------------------------------------------------
% OUTPUT
%   scores (numeric): Computed scores for the model.
%--------------------------------------------------------------------------
% EXAMPLES
%  See: test_compute_scores.m
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

if ~exist('score_type', 'var')
    score_type = 'effective';
end

y = y';
nsubj = size(y, 1);

% Obtain the derivative and variance evaluated at the fitted values
[Dhat0, Vhat0] = get_par_expo_fam(fitted_values_0, family, link);

% Initialize the weights matrix to be ones for now - this may change!
W = (Dhat0 .^ 2) ./ Vhat0;
sqrtW = sqrt(W);

% Compute the inverse square root of V
sqrtinvVvect = Vhat0 .^ (-0.5);

% Obtain the residuals
null_residuals = y - fitted_values_0;

sqrtinvVvect_times_residuals = sqrtinvVvect .* null_residuals;

if strcmp(score_type, 'effective')
    A = Z'.*sqrtW';  % Calculate Z transpose times diag(sqrtW)
    XTsqrtW = X'.*sqrtW';
    H = A' * inv(A * A') * A;
    scores = (XTsqrtW * (eye(nsubj) - H)).* (sqrtinvVvect_times_residuals' / (nsubj^0.5));
end

end
