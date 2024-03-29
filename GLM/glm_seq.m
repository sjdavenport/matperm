function [betahat, fitted_values, perfect_separations, pvalues, std_errors] = ...
                          glm_seq(y, X, family, linkfn, progress_message)
% glm_seq runs the generalized linear model sequentially at each voxel in
% the array y.
% Inputs:
%   y: a matrix of size (n_voxels, n_samples) representing the dependent variable.
%   X: a matrix of size (n_samples, n_parameters) representing the independent variables.
%   distbn: a string specifying the distribution of the response variable,
%       one of {'Binomial', 'Gamma', 'Gaussian', 'Inverse Gaussian', 'Negative Binomial', 'Poisson'}.
%   linkfn: a string specifying the link function, one of {'log', 'logit', 'probit', 'cauchy', 'cloglog', 'identity', 'inverse'}.
%   progress_message: a format string to display progress, e.g., 'Progress: %d / %d\n'.
% Outputs:
%   betahat: a matrix of size (n_parameters, n_voxels) containing the estimated coefficients for each voxel.
%   fitted_values: a matrix of size (n_samples, n_voxels) containing the fitted values.
%   perfect_separations: a vector indicating voxels with perfect separation.
%   pvalues: a matrix of size (n_parameters, n_voxels) containing p-values for coefficients.
%   std_errors: a matrix of size (n_parameters, n_voxels) containing standard errors for coefficients.
%
% Examples:
%  See test_glm_seq.m
% Note: This function performs a generalized linear model (GLM) analysis 
% for each voxel in the input data.

if ~exist('progress_message', 'var')
   progress_message = 'GLM progress: '; 
end

y = y'; % Transpose y

% Calculate the number of voxels
nvox = size(y, 2);

% Initialize arrays to store results
s_X = size(X, 2);
betahat = zeros(nvox, s_X);
pvalues = zeros(nvox, s_X);
std_errors = zeros(nvox, s_X);
llf = zeros(nvox, 1);
fitted_values = zeros(nvox, size(y, 1));
conf_ints = zeros(nvox, s_X, 2);
perfect_separations = zeros(nvox, 1);

% warning('off', 'MATLAB:log:logOfZero');

warning("")
% Loop over the columns of y
for i = 1:nvox
    % Measure progress
    if ~isempty(progress_message)
        loader(i, nvox, progress_message);
    end
    % Fit a GLM for the i-th column of y
    try
        [model, ~, s] = glmfit(X, y(:, i), family, 'link', linkfn, 'constant', 'off');
        betahat(i,:) = model';
        fitted_values(i,:) = glmval(model, X, linkfn, 'constant', 'off')';
        pvalues(i,:) = s.p';
        std_errors(i,:) = s.se';
%         llf(i) = model(1); % Log-likelihood
%         conf_ints(:, i, :) = coefCI(model);
    catch
        perfect_separations(i) = 1;
    end
    if strcmp('Iteration limit reached.', lastwarn)
        perfect_separations(i) = 1;
        warning("")
    end
end

% warning('on', 'MATLAB:log:logOfZero');

perfect_separations_second_loc = find(sum(fitted_values == 0) > 0);
perfect_separations(perfect_separations_second_loc) = 1;

perfect_separations_third_loc = find(sum(fitted_values == 1) > 0);
perfect_separations(perfect_separations_third_loc) = 1;

% Return the results
end
