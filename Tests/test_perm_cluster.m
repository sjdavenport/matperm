%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the perm_cluster function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example


%% %% 2D Examples
%% Simple 2D example - cluster size vs tfce
connectivity_criterion = 8; H = 2; E = 0.5;
dim = [50,50]; nsubj = 50; FWHM = 0;
Sig = 0.5*peakgen(1, 10, 8, dim);
% Sig = zeros(dim); Sig(25:26,25) = 3;
% Sig = 0.5*square_signal(dim, 4, {[25,20], [25,30]} );
data = wfield(dim, nsubj).field + Sig;
threshold_tfce = perm_tfce(data, ones(dim), H, E, connectivity_criterion);
tstat_orig = mvtstat(data);
tfce_tstat = tfce(tstat_orig, H, E, connectivity_criterion);
CDT = 2.3;
threshold_cluster = perm_cluster(data, ones(dim), CDT, connectivity_criterion);

[number_of_clusters, occurences, sizes, index_locations] = numOfConComps(tstat_orig, CDT, connectivity_criterion);
surviving_cluster_im = cluster_im( dim, index_locations, threshold_cluster );

%%
subplot(2,2,1)
surf(tstat_orig)
title('Original t-stat')
view([-3, 15.58])
subplot(2,2,2)
surf(tfce_tstat)
view([-3, 15.58])
title('TFCE statistic')
subplot(2,2,3)
imagesc(surviving_cluster_im)
title('Cluster extext inference: CDT = 2.3')
subplot(2,2,4)
imagesc(tfce_tstat > threshold_tfce)
title('TFCE: H = 2, E = 0.5')
% axis square
fullscreen
export_fig('C:\Users\12SDa\global\TomsMiniProject\Latex\MyPapers\TFCE_vs_ClusterExtent\round_sig.png')

%% Simple 2D example - cluster size vs copesets
connectivity_criterion = 8; H = 2; E = 0.5;
dim = [50,50]; nsubj = 50; FWHM = 3;
Sig = 0.5*peakgen(1, 10, 8, dim);
% Sig = zeros(dim); Sig(25:26,25) = 3;
% Sig = 0.5*square_signal(dim, 4, {[25,20], [25,30]} );
data = wfield(dim, nsubj).field + Sig;
data = fconv(data, FWHM, 2);

% Run clustersize inference
CDT = 2.3;
threshold_cluster = perm_cluster(data, ones(dim), CDT, connectivity_criterion);

tstat_orig = mvtstat(data);
[number_of_clusters, occurences, sizes, index_locations] = numOfConComps(tstat_orig, CDT, connectivity_criterion);
surviving_cluster_im = cluster_im( dim, index_locations, threshold_cluster );

% Run SSS cope set inference
nboot = 1000;
c = 0;
[lower_set, upper_set, std_multipler] = sss_cope_sets(data, ones(dim), c, nboot);

%%
subplot(1,3,1)
surf(tstat_orig)
title('Original t-stat')
axis square
view([-3, 15.58])
subplot(1,3,2)
imagesc(surviving_cluster_im)
title('Cluster extext inference: CDT = 2.3')
axis square
subplot(1,3,3)
cope_display( lower_set, upper_set, mean(data,3), c);
title('Cope set image')
% axis square
fullscreen
% export_fig('C:\Users\12SDa\global\TomsMiniProject\Latex\MyPapers\TFCE_vs_ClusterExtent\round_sig.png')

%% Comparing cluster size inference and scopes
[ lower_band, upper_band ] = scopes( data, 1000, 0.05 );
c_vec = 0:0.25:0.5;
subplot(1,length(c_vec)+1,1)
imagesc(surviving_cluster_im)
axis square
title('Cluster size inference')
for I = 1:length(c_vec)
    subplot(1, length(c_vec)+1,I+1)
    cope_display( upper_band > c_vec(I), lower_band > c_vec(I), mean(data,3), c_vec(I));
    title(['SCOPES, c = ', num2str(c_vec(I))])
end
%% SCOPES
[ lower_band, upper_band ] = scopes( data, 1000, 0.05 );
subplot(2,1,1)
cope_display( upper_band > c, lower_band > c, mean(data,3), c);
subplot(2,1,2)
cope_display( lower_set, upper_set, mean(data,3), c);

%%
imagesc(upper_band)
%% Nearby square signals
dim = [50,50]; nsubj = 50; FWHM = 0;
Sig = 0.25*peakgen(1, 10, 8, dim);
Sig = 0.5*square_signal(dim, 4, {[25,20], [25,30]} );
data = wfield(dim, nsubj);
data.field = data.field + Sig;
tstat = convfield_t(data, FWHM);
tstat_tfce = tfce(tstat.field,2,0.5,8,0.05);

subplot(1,2,1)
surf(tstat.field)
title('Original tstat')
view([-14,15])
subplot(1,2,2)
surf(tstat_tfce)
title('TFCE')
view([-14,15])
fullscreen

%%
CDT = 1;
threshold_cluster = perm_cluster(data, ones(dim), CDT, connectivity_criterion);

[number_of_clusters, occurences, sizes, index_locations] = numOfConComps(tstat_orig, CDT, connectivity_criterion);
surviving_cluster_im = cluster_im( dim, index_locations, threshold_cluster );


%% %% 3D Examples
%% Simple 3D example