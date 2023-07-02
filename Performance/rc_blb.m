function fpr = rc_blb(  spfn, nsubj, subset_size, n_subsamples, n_mc, niters, alpha )
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
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'alpha', 'var' )
   % Default value
   alpha = 0.05;
end

if ~exist( 'niters', 'var' )
   % Default value
   niters = 1000;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
fpr = 0;
for b = 1:niters
    modul(b,100)
    data = spfn(nsubj);
    [ bootstrap_samples ] = blb( data, subset_size, n_subsamples, n_mc );
    threshold = prctile(bootstrap_samples, 100*(1-alpha) );
    if bootstrap_samples(1) > threshold
        fpr = fpr + 1;
    end
end
fpr = fpr/niters;

end

