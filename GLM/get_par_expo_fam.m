function [Dhat, Vhat] = get_par_expo_fam(fitted_values, family, link)
% GET_PAR_EXPO_FAM - Compute the derivative and variance for exponential family.
%
% Parameters:
%   fitted_values (array-like): The fitted values for which to compute the
%                               derivative and variance.
%   family (char): The exponential family to use, can be 'gaussian',
%                  'binomial', 'poisson', or 'gamma'.
%   link (char): The link function to use, can be 'identity', 'log', 'logit',
%                'cloglog', 'probit', or 'cauchit'.
%
% Returns:
%   Dhat (numeric): The derivative of the link function at the fitted values.
%   Vhat (numeric): The variance function at the fitted values.
%
% Examples:
%   [Dhat, Vhat] = get_par_expo_fam(fitted_values, 'binomial', 'logit');
%   % Compute Dhat and Vhat for a binomial distribution with logit link.
%
% See also: get_linkfn (if applicable).

fitted_values = double(fitted_values);
if strcmp(family, 'Gaussian')
    var_fn = @(x) 0 * x + 1;
elseif strcmp(family, 'Binomial') || strcmp(family, 'binomial')
    var_fn = @(x) x .* (1 - x);
elseif strcmp(family, 'Poisson')
    var_fn = @(x) x;
elseif strcmp(family, 'Gamma')
    var_fn = @(x) x .* x;
else
    error('The family must be one of the specified ones');
end

% Compute the variance function at the fitted values
Vhat = var_fn(fitted_values);

% Get the derivative of the link function
[link_fn, ~, link_fn_inv_deriv] = get_linkfn(link);

% Compute the derivative of the link at the fitted values
Dhat = link_fn_inv_deriv(link_fn(fitted_values));
end
