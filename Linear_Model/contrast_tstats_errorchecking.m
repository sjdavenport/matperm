function [ contrast_matrix, nsubj, n_params ] = contrast_tstats_errorchecking( ...
                                 data, design_matrix, contrast_matrix )
% CONTRAST_TSTATS_ERRORCHECKING( data, design, contrast_matrix )
% performs error checking on the contrast data to ensure that it has the
% right dimensions.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% data:  an N by V matrix (N is the number of subjects and V the number of
%           voxels
% design:    an N by p matrix giving the covariates (p being the number of
%            parameters
% contrast_matrix: an L by p matrix corresponding to the contrast matrix, 
%                  such that which each row is a contrast vector (where L 
%                  is the number of constrasts)
%--------------------------------------------------------------------------
% OUTPUT
% contrast_matrix: an L by p matrix corresponding to the contrast matrix, 
%                  such that which each row is a contrast vector (where L 
%                  is the number of constrasts)
% nsubj:    an integer giving the number of subjects
% n_params: an integer giving the number of parameters in the model
%--------------------------------------------------------------------------
% EXAMPLES
% Dim = [3,3]; N = 30; categ = zeros(1, N);
% X = group_design(categ); C = 1; lat_data = wfield(Dim,N);
% [ contrast_matrix, nsubj, n_params ] = contrast_tstats_noerrorchecking(lat_data, X, C)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% s_contrast_matrix = size(contrast_matrix);

%%  Main Function Loop
%--------------------------------------------------------------------------
%% Error Checking
% Calculate the number of parameters in C
n_contrast_params = size(contrast_matrix,2); % parameters

% Calculate the number of parameters p and subjects N
nsubj = size(design_matrix,1); % subjects
n_params = size(design_matrix,2); % parameters

% Ensure that the dimensions of X and C are compatible
if n_params ~= n_contrast_params
    error('The dimensions of design matrix and contrast_matrix do not match')
end

s_lat_data = size(data);
% Ensure that the dimensions of X and lat_data are compatible
if nsubj ~= s_lat_data(end)
    error('The number of subjects in design matrix and lat_data do not match')
end

end

