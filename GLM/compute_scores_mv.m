function [scores, pslocs] = compute_scores_mv(y, Z, X, family, link, score_type, progress_message)
% COMPUTE_SCORES_MV(y, Z, X, fitted_values_0, family, link, score_type)
% computes the score contributions for a generalized linear model.
%--------------------------------------------------------------------------
% INPUT
%   y (numeric): A matrix of shape (nvoxels, nsubj) representing the dependent variable.
%   Z (numeric): An array of shape (nsubj, n_nuisance_parameters) representing the nuisance variables.
%   X (numeric): An array of shape (nsubj, n_parameters) representing the independent variables.
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
% Copyright (C) - 2023 - Samuel Davenport
%--------------------------------------------------------------------------

%%  Set optional variables
%--------------------------------------------------------------------------
if ~exist('score_type', 'var')
    score_type = 'effective';
end

if ~exist('progress_message', 'var')
    progress_message = 'Fitting null glm model, progress:';
end

%% Get the important constants
%--------------------------------------------------------------------------
nvox = size(y, 1);
nsubj = size(y, 2);

if size(X,1) ~= nsubj
    error('X must have nsubj rows');
end
if size(Z,1) ~= nsubj
    error('Z must have nsubj rows');
end

p = size(X, 2);

%%  Main Function Loop
%--------------------------------------------------------------------------
% Compute null model
if strcmp(score_type, 'firth')
    [~, fitted_values_0, psZ] = firth_regression_seq(y,Z,progress_message);
    score_type = 'effective';
else
    [~, fitted_values_0, psZ] = ...
        glm_seq(y, Z, family, link, progress_message);
end

% Combine the perfect separations locations into a single vector
pslocs = find(psZ);
% Compute the scores
if ~isempty(progress_message)
    disp('Computing scores...');
end
scores = zeros(nvox, p, nsubj);
for I = setdiff(1:nvox, pslocs)
    yvox = y(I,:);
    fitted_0_vox = fitted_values_0(I,:);
    scores(I, :, :) = compute_scores(yvox, Z, X, fitted_0_vox', family, link, score_type);
end

end
