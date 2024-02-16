function [ tstat_array, residuals, Cbetahat, betahat, sigmahat ] = ...
             contrast_tstats( data, design, contrast_matrix, check_error )
% CONTRAST_TSTATS( lat_data, design, contrast_matrix, check_error )
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% 	data:  an array of size (Dim, N) or an object of class field
%      	giving the data where Dim is the spatial dimension and N is the
%       number of subjects if a field then the fibersize must be 1 and 
%       the final dimension must be the number of subjects
% 	design: a matrix of size (N,p)
%       giving the covariates (p being the number of parameters)
% 	contrast_matrix: a matrix of size (L,p)
%       corresponding to the contrast matrix, such that which each row is a
%       contrast vector (where L is the number of constrasts)s
% Optional
%   check_error:  0/1, determining whether to perform error checking or not
%       (not always necessary e.g. during a permutation loop etc) default
%       is 1 i.e. to perform error checking
%--------------------------------------------------------------------------
% OUTPUT
%   tstat_field: can object of class field which has spatial size the same
%                as input data and fibersize equal to the number of contrasts
%   residuals: a matrix of size (dim, nsubj) containing the residuals i.e. 
%               residuals(..., I) provides the imagewise residuals for the 
%               Ith subject
%--------------------------------------------------------------------------
% EXAMPLES
% % One-sample t-statistic
% Dim = [3,3]; N = 30; categ = zeros(1, N);
% X = group_design(categ); C = 1; lat_data = wfield(Dim,N);
% tstat_array = contrast_tstats(lat_data.field, X, C)
% %Compare to mvtstat:
% tstat_array
% mvtstat(lat_data.field)
%
% % Several contrasts
% age = randi([40,60], nsubj, 1);
% sex = rand(nsubj, 1) > 0.6;
% Z = [ones(nsubj,1), sex, age];
% X = randn(nsubj,1);
% design_matrix = [X,Z];
% nvox = 1000;
% X_coeff = 2;
% data = 5 + 0.12*(age'-50) + 0.2*sex' + X_coeff*X' + randn(nvox, nsubj);
% tstats = contrast_tstats( data, design_matrix, [1,0,0,0] )
% histogram(tstats)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'check_error', 'var' )
    % Default value
    check_error = 1;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% Error check the inputs
if check_error == 1
    contrast_matrix = contrast_tstats_errorchecking(data, design, contrast_matrix);
end

% Having now run error checking calculate the contrast t-statistics
[ tstat_array, residuals, Cbetahat, betahat, sigmahat ] = contrast_tstats_noerrorchecking(data, design, contrast_matrix);

end

