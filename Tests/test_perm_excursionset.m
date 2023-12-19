%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the perm_excursionset function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2D example
dim = [50,50]; nsubj = 50; 
% Sig = 0.5*peakgen(1, 10, 8, dim);
% Sig = zeros(dim); Sig(25:26,25) = 3;
% Sig = 0.35*square_signal(dim, 4, {[25,20], [25,30]} );
Sig = 0.35*square_signal(dim, 10, {[25,14], [25,36]} );
data = wfield(dim, nsubj).field + Sig;
% FWHM = 0; 
FWHM = 5;
data = fconv(data, FWHM, 2);
tstat_orig = mvtstat(data);
CDT = 2.3;
CDT = 3.1;
threshold_excursionset = perm_excursionset(data, ones(dim), CDT);
threshold_cluster = perm_cluster(data, ones(dim), CDT);

%%
index_locations = find(tstat_orig > CDT);
if length(index_locations) > threshold_excursionset
    imagesc(tstat_orig > CDT)
end

%%


%% %% 2D Examples
%% Simple 2D example


%% %% 3D Examples
%% Simple 3D example