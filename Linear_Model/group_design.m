function [ design ] = group_design( categ )
% A function to compute the covariate matrix for a given set of categories
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% categ:  a tuple of integers of length N
%         where N is the number of subjects). Each entry is number of the category
%         that a given subject belongs to (enumerated from 0 to ncateg - 1)
%         E.g: (0,1,1,0) corresponds to 4 subjects, 2 categories and
%                  (0,1,2,3,3,2) corresponds to 6 subjects and 4 categories
%--------------------------------------------------------------------------
% OUTPUT
% design: a design matrix that can be used to assign the correct categories
%--------------------------------------------------------------------------
% EXAMPLES
% categ = [0,1,1,0]; group_design(categ)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % Default value
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% Calculate the number of subjects
nsubj = length(categ);

% Calculate the number of parameters i.e. the number of distinct groups
n_params = length(unique(categ));

% Ensure that the number of categories is not too high!
if max(categ) > n_params - 1
    error("the maximum category number should not exceed one minus the number of categories")
end

% Initialize the design matrix
design = zeros(nsubj,n_params);

% Set the elements of the design matrix by assigning each subject a category
for I = 1:nsubj
	design(I, categ(I) + 1) = 1;
end

end

