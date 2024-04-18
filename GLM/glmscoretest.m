function [ out ] = glmscoretest(y, Z, X, family, link)
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
if ~exist( 'opt1', 'var' )
   % Default value
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
[~, fitted_values_0, psZ] = glm_seq(y, Z, family, link, []);
[scores, XTsqrtWIminusH, sqrtinvVvect_times_residuals] = ...
        compute_scores(y, Z, X, fitted_values_0, family, link, 'effective')

score_test = mean(scores);

end

