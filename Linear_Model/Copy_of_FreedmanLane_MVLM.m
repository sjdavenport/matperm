function [ tstat_perm ] = FreedmanLane( Y, X, Z, contrast_vector, niters )
% NEWFUN serves as a function template.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   Y
%   X
%   Z
%   contrast_vector
%   niters
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Calculate the number of subjects
nsubj = size(X,1);

% Calculate H_Z = Z(Z^TZ)^(-1)Z^T
H_Z = Z*pinv(Z);
% Calculate R_Z = I - H_Z
R_Z = eye(nsubj) - H_Z;

% Initialize a vector to store the permuted t-statistics
tstat_perm = zeros(1, niters);

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % Default value
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------

% Fit the linear model on the original data
fit_mvlm_real = MVlm_multivar( [X,Z], Y, contrast_vector, 1 );

% Set the first permuted t-stat to be the one from the fit
tstat_perm(1) = fit_mvlm_real.tstat;

%% Permutation loop
for I = 2:niters
    % Generate a random sample
    rand_perm = randsample(1:10, 10, 0);
    
    % Shuffle the rows of the R_Z matrix
    R_Z_perm = R_Z(rand_perm,:);
    Yperm = (R_Z_perm + H_Z)*Y;
    fit_mvlm_perm = MVlm_multivar( [X,Z], Yperm, contrast_vector, 1 );
    tstat_perm(I) = fit_mvlm_perm.tstat;
end

end

