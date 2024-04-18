function [betahat, fitted_values, noconvergence_locs] = firth_regression_seq( y, X, progress_message)
% NEWFUN
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% Copyright (C) - 2023 - Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist('progress_message', 'var')
   progress_message = 'Firth regression progress: '; 
end

y = y'; % Transpose y

% Calculate the number of voxels
nvox = size(y, 2);

s_X = size(X, 2);

betahat = zeros(nvox, s_X);
fitted_values = zeros(nvox, size(y, 1));
noconvergence_locs = zeros(nvox, 1);
%%  Main Function Loop
%--------------------------------------------------------------------------
for i = 1:nvox
    % Measure progress
    if ~isempty(progress_message)
        loader(i, nvox, progress_message);
    end
    
    % Fit a GLM for the i-th column of y
    
    try
    [betahat(i,:) , fitted_values(i,:) ] = firth_regression(y(:, i), X);
    catch
        noconvergence_locs(i) = 1;
    end
end

end

