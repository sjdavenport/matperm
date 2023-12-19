function [scores, pslocs] = compute_scores_mv(y, Z, X, distbn, link)
    % Compute null model
    [~, fitted_values_0, psZ] = ...
                          glm_seq(y, X, distbn, link, 'Fitting null glm model, progress:');

    % Combine the perfect separations locations into a single vector
    pslocs = find(psZ);

    % Get the important constants
    nvox = size(y, 1);
    nsubj = size(y, 2);
    p = size(X, 2);

    % Compute the scores
    disp('Computing scores...');
    scores = zeros(nvox, p, nsubj);
    for I = setdiff(1:nvox, pslocs)
        yvox = y(I,:);
        fitted_0_vox = fitted_values_0(:, i);
        scores(i, :, :) = compute_scores(yvox, Z, X, fitted_values_0, family, link, 'effective');
    end
end
