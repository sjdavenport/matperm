function fp_power( effect_size, nsubj, domeanperm, numberofblocks, nperm, nvox, niters, demean )
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
if ~exist( 'domeanperm', 'var' )
   % Default value
   domeanperm = 1;
end

if ~exist( 'niters', 'var' )
   % Default value
   niters = 1000;
end 

if ~exist( 'demean', 'var' )
   % Default value
   demean = 0;
end 

if ~exist( 'nperm', 'var' )
   % Default value
   niters = 1000;
end 

if ~exist('numberofblocks', 'var')
    numberofblocks = [5,10,20,50,100,250,1000];
end


%%  Main Function Loop
%--------------------------------------------------------------------------
fpr_fp = cell(1,2);

% Initialize the storage vectors
for init_randomize = 0:1
    fpr_fp{init_randomize+1} = zeros(1, length(numberofblocks));
end

for I = 1:length(numberofblocks)
    I
    spfn = @(n) randn(n, nvox)' + effect_size;
    for init_randomize = 0:1
        for demean = 0:1
            fpr_fp{init_randomize+1, demean+1}(I) = rc_fastperm( spfn, nsubj, numberofblocks(I), nperm, domeanperm, init_randomize, demean, niters);
        end
    end
end

mploc = 'C:\Users\12SDa\davenpor\davenpor\Toolboxes\matperm\';
if effect_size == 0
    if domeanperm == 1
        save([mploc,'Performance/FPR/fpr_vs_nblocks_mean_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox)], 'fpr_fp', 'niters', 'numberofblocks')
    else
        save([mploc,'Performance/FPR/fpr_vs_nblocks_tstat_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox)], 'fpr_fp', 'niters', 'numberofblocks')
    end
else
    if domeanperm == 1
        save([mploc,'Performance/Power/power_vs_nblocks_mean_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox)], 'fpr_fp', 'niters', 'numberofblocks')
    else
        save([mploc,'Performance/Power/power_vs_nblocks_tstat_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox)], 'fpr_fp', 'niters', 'numberofblocks')
    end
end

end

