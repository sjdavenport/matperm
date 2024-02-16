function [beta, fitted_values, bse, fitll] = firth_regression(y, X, start_vec, step_limit, convergence_limit)
    if nargin < 3
        start_vec = zeros(size(X, 2), 1);
    end
    if nargin < 4
        step_limit = 1000;
    end
    if nargin < 5
        convergence_limit = 0.0001;
    end
    % 
    % % Convert y to categorical variable
    % y_cat = double(categorical(y, [0, 1]));

    %b = glmfit(X, y, 'binomial', 'link', 'logit', 'constant', 'off');

    beta_iterations = zeros(step_limit, size(X, 2));
    beta_iterations(1, :) = start_vec';

    for i = 1:step_limit
        pi = glmval(beta_iterations(i, :)', X, 'logit','constant', 'off');
        % pi = glmval(beta_iterations(i, :)', X, 'identity','constant', 'off');
        W = diag(pi .* (1 - pi));
        var_covar_mat = pinv(X' * W * X);

        % build hat matrix
        rootW = sqrt(W);
        H = X' * rootW;
        H = var_covar_mat * H;
        H = rootW * X * H;

        % penalised score
        U = X' * (y - pi + diag(H) .* (0.5 - pi));
        new_beta = beta_iterations(i, :)' + var_covar_mat * U;

        % step halving
        j = 0;
        while firth_likelihood(new_beta,X,y) > firth_likelihood(beta_iterations(i, :)', X, y)
            new_beta = beta_iterations(i, :)' + 0.5 * (new_beta - beta_iterations(i, :)');
            j = j + 1;
            if j > step_limit
                fprintf('Firth regression failed. Try increasing step limit.\n');
                return;
            end
        end

        beta_iterations(i + 1, :) = new_beta';
        if i > 1 && norm(beta_iterations(i, :)' - beta_iterations(i - 1, :)') < convergence_limit
            break;
        end
    end

    return_fit = [];

    if norm(beta_iterations(i, :)' - beta_iterations(i - 1, :)') >= convergence_limit
        fprintf('Firth regression failed to converge.\n');
    else
        % Calculate stats
        fitll = -firth_likelihood(beta_iterations(i, :)', X, y);
        beta = beta_iterations(i, :)';
        bse = sqrt(diag(pinv(X' * W * X)));
    end

    fitted_values = sigmoid(X*beta);
end

