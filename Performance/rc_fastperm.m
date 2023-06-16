function [ fpr ] = rc_fastperm( spfn, nsubj, nblocks, nperm, domeanperm, init_randomize, demean, niters, alpha )
% RF_FASTPERM tests the fpr for the fast permutation method
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  spfn:     a function handle that generates data such that spfn(nsubj) is 
%            an nvox by nsubj array
%  nsubj:    an integer giving the number of subjects
%  nblocks:  the number of blocks with which to divide the data for the fast
%            permutation    
% Optional
%  nperm     an integer giving the number of permutations
%  init_randomize  0/1 whether or not to randomly sign-flip initially, 
%            default is 1 i.e. to do an initial sign-flip
%  niters    an integer giving the number of iterations to do, default is 
%            1000
%  alpha     a number between 0 and 1 giving the alpha level at which to
%            control the fpr, default is 0.05
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

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'nperm', 'var' )
   % Default value
   nperm = 1000;
end

if ~exist( 'domean', 'var' )
   % Default value
  	domeanperm = 0;
end

if ~exist( 'alpha', 'var' )
   % Default value
   alpha = 0.05;
end

if ~exist( 'niters', 'var' )
   % Default value
   niters = 1000;
end

if ~exist( 'demean', 'var' )
   % Default value
   demean = 0;
end

if ~exist( 'init_randomize', 'var' )
   % Default value
   init_randomize = 1;
end
%%  Main Function Loop
%--------------------------------------------------------------------------
fpr = 0;
for b = 1:niters
    modul(b,100)
    data = spfn(nsubj);
    if domeanperm
        [threshold, vec_of_maxima, ~ ] = fastperm_mean( data, nblocks, alpha, nperm, demean, 0, 0, init_randomize);
    else
        [threshold, vec_of_maxima, ~ ] = fastperm( data, nblocks, alpha, nperm, 0, 0, init_randomize);
    end
    
    if vec_of_maxima(1) > threshold
        fpr = fpr + 1;
    end
end
fpr = fpr/niters;

end

