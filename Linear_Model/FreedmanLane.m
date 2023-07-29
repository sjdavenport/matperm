function [ vec_of_maxima ] = FreedmanLane( Y, X, Z, contrast_vector, nperm )
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
% nsubj = 100;
% Y = randn(2, nsubj);
% Z = randn(nsubj,1);
% X = randn(nsubj,1);
% vec_of_maxima = FreedmanLane( Y, X, Z, [1,0], 100 );
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'nperm', 'var' )
   % Default value
   nperm = 1000;
end


%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Calculate the number of subjects
nsubj = size(X,1);

% Calculate H_Z = Z(Z^TZ)^(-1)Z^T
H_Z = Z*pinv(Z);
% Calculate R_Z = I - H_Z
R_Z = eye(nsubj) - H_Z;

% Initialize a vector to store the permuted t-statistics
vec_of_maxima = zeros(1, nperm);

%%  Main Function Loop
%--------------------------------------------------------------------------

% Fit the linear model on the original data
tstat_orig = contrast_tstats( Y, [X,Z], contrast_vector );
% tstat_orig = contrast_tstats( Y*R_Z', [X,Z], contrast_vector ); this is
% the same as the above interestingly!!!

% Set the first permuted t-stat to be the one from the fit
vec_of_maxima(1) = max(tstat_orig);

random_berns = 2*(binornd(1,0.5,nsubj,nperm )-1/2);

% HZY = H_Z*Y';

%% Permutation loop
for I = 1:nperm-1
    % Generate a random sample
%     rand_perm = randsample(nsubj, nsubj, 0);
%     
%     % Shuffle the rows of the R_Z matrix
%     R_Z_perm = R_Z(rand_perm,:);
    R_Z_perm = R_Z;
    
    random_berns_for_iter = random_berns(:, I);
    random_sample_negative = find(random_berns_for_iter < 0);
    R_Z_perm(random_sample_negative, :) = -R_Z_perm(random_sample_negative, :);
    
%     Yperm = (R_Z_perm*Y' + HZY)'; % could be made more efficient here!!
    Yperm = Y*R_Z_perm';
%     fit_mvlm_perm = MVlm_multivar( [X,Z], Yperm, contrast_vector, 1 );
    tstat_perm_I = contrast_tstats( Yperm, [X,Z], contrast_vector );

    vec_of_maxima(I+1) = max(tstat_perm_I);
end

end

