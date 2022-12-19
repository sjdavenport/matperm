function [ tstat_array, residuals, Cbetahat, betahat, sigmahat ] = ...
      contrast_tstats_noerrorchecking(data, design_matrix, contrast_matrix)
% A function to compute the voxelwise t-statistics for a set of contrasts
% but with no error checking. This is design_matrixed for input into permutation 
% so you do not have to run the error checking every time.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% data:  an object of class field giving the data for N subjects on 
%        which to calculate the contrasts
% design_matrix:    an N by p matrix giving the covariates (p being the 
%                   number of parameters
% contrast_matrix: an L by p matrix corresponding to the contrast matrix, 
%                  such that which each row is a contrast vector (where L 
%                  is the number of constrasts)
%--------------------------------------------------------------------------
% OUTPUT
% tstat_field: can object of class field which has spatial size the same
%                as input data and fibersize equal to the number of contrasts
% residuals: a matrix of size (dim, nsubj) containing the residuals i.e. 
%            residuals(..., I) provides the imagewise residuals for the 
%            Ith subject
% Cbeta_field: an object of class field which has shape (dim, ncontrasts), 
%            for each contrast this gives the c^Tbetahat image.
%--------------------------------------------------------------------------
% EXAMPLES
% % One-sample t-statistic
% Dim = [3,3]; N = 30; categ = zeros(1, N);
% X = group_design_matrix(categ); C = 1; lat_data = wfield(Dim,N);
% tstat = contrast_tstats_noerrorchecking(lat_data.field, X, C)
% %Compare to mvtstat:
% tstat
% mvtstat(lat_data.field)
% %--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Main Function Loop
%--------------------------------------------------------------------------
n_contrasts = size(contrast_matrix,1);  % number of constrasts
D = length(size(data)) - 1;

% Calculate the number of parameters p and subjects N
nsubj = size(design_matrix,1); % subjects
n_params = size(design_matrix,2); % parameters

xtx_inv = inv(design_matrix'*design_matrix);

% Calculate betahat everywhere
betahat = marray(xtx_inv * design_matrix', data);
% betahat = xtx_inv * design_matrix' * lat_data;

% Calculate the residual forming matrix
rfmate = eye(nsubj) - design_matrix * xtx_inv * design_matrix';

% Compute the estimate of the variance via the residuals (I-P)Y
% residuals = rfmate * lat_data;
residuals = marray(rfmate, data);

%Square and sum over subjects to calculate the variance
% This assumes that X has rank p!
% sigmahat = (sum(residuals.field.^2,lat_data.D + 1)*(1/(nsubj-n_params))).^(1/2);
sigmahat = (sum(residuals.^2,D + 1)*(1/(nsubj-n_params))).^(1/2);

% Compute the t-statistics
if isequal(size(contrast_matrix), [1,1])
    Cbetahat = contrast_matrix*betahat;
else
    Cbetahat = marray(contrast_matrix, betahat);
end
tstat_array = Cbetahat./sigmahat;

% Create a variable index
variable_index = repmat( {':'}, 1, D );

% Scale by the scaling constants to ensure variance 1
for l = 1:n_contrasts
    scaling_constant = (contrast_matrix(l,:) * xtx_inv * contrast_matrix(l,:)').^(1/2);
	tstat_array(variable_index{:},l) = tstat_array(variable_index{:},l)/scaling_constant;
end

% tstat_field = Field(tstats, lat_data.mask);

end

