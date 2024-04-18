function [likelihood, log_likelihood, hessian] = firth_likelihood(beta, X, y)
    % y = y';

    % Compute predicted probabilities
    p = 1 ./ (1 + exp(-X * beta));
    
    % Compute log likelihood
    log_likelihood = sum(y .* log(p) + (1 - y) .* log(1 - p));
    
    % Compute Hessian matrix
    n = length(y);
    hessian = zeros(length(beta), length(beta));
    
    for k = 1:n
        hessian = hessian - p(k) * (1 - p(k)) * X(k, :)' * X(k, :);
    end

    likelihood = -(log_likelihood + 0.5 * log(det(-hessian)));
end