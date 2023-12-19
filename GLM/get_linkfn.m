function [link_fn, link_fn_inv, link_fn_inv_deriv] = get_linkfn(link)
    % GET_LINKFN - Returns link functions for the specified link type.
    %
    % Syntax:
    %   [link_fn, link_fn_inv, link_fn_inv_deriv] = get_linkfn(link)
    %
    % Inputs:
    %   link - A string specifying the link function type. Choose from:
    %          'identity', 'log', 'logit', 'cloglog', 'probit', 'cauchit'.
    %
    % Outputs:
    %   link_fn - A function handle for the specified link function.
    %   link_fn_inv - A function handle for the inverse of the link function.
    %   link_fn_inv_deriv - A function handle for the derivative of the
    %                      inverse of the link function.
    %
    % Example:
    %   [lf, lfi, lfiv] = get_linkfn('logit');
    %   y = lf(0.5)           % Compute the link function value
    %   x = lfi(0.5)         % Compute the inverse link function value
    %   dx = lfiv(0.5)       % Compute the derivative of the inverse link function
    
    if strcmp(link, 'identity')
        link_fn_inv = @(x) x;
        link_fn_inv_deriv = @(x) 1;
        link_fn = @(x) x;
    elseif strcmp(link, 'log')
        link_fn_inv = @(x) exp(x);
        link_fn_inv_deriv = @(x) exp(x);
        link_fn = @(x) log(x);
    elseif strcmp(link, 'logit')
        link_fn_inv = @(x) exp(x) ./ (1 + exp(x));
        link_fn_inv_deriv = @(x) exp(x) ./ (1 + exp(x)).^2;
        link_fn = @(x) log(x ./ (1 - x));
    elseif strcmp(link, 'cloglog')
        link_fn_inv = @(x) 1 - exp(-exp(x));
        % You seem to be missing 'link_fn_deriv' in your original Python code.
    elseif strcmp(link, 'probit')
        link_fn_inv = @(x) exp(-x.^2 / 2) / sqrt(2 * pi);
        % You seem to be missing 'link_fn_deriv' in your original Python code.
    elseif strcmp(link, 'cauchit')
        link_fn_inv = @(x) 1 / (pi * (1 + x.^2));
        % You seem to be missing 'link_fn_deriv' in your original Python code.
    else
        error('The link function must be one of the specified ones');
    end
end
