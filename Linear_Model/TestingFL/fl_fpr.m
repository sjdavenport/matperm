nsubj = 10;
nvox = 100;
% Y = randn(nvox, nsubj);
% Z = randn(nsubj,1);
% % Z 
% X = randn(nsubj,1);
% % Y = Z' + X' + randn(nvox,nsubj);
% Y = Z' + randn(nvox,nsubj);

niters = 1000;
fpr = 0;
for I = 1:niters
    I
    modul(I,50)
    Z = randn(nsubj,1);
    X = randn(nsubj,1);
%     Y = Z' + X' + randn(nvox,nsubj);
%     Y = Z' + rlap( 3, [nvox,nsubj] );
    sample = randsample(198,nsubj,0);
    vox_sample = randsample(size(RSD_data_vec, 1), nvox, 0);
    Y = 100*Z' + -RSD_data_vec(vox_sample, sample);
%     Y = Z' + X' + RSD_data_vec(vox_sample, sample);
vec_of_maxima = FreedmanLane( Y, X, Z, [1,0], 1000 );
% vec_of_maxima = FreedmanLane_Tom(Y',X,Z, 1000);
alpha_level = 0.05;
thresh = prctile(vec_of_maxima, 100*(1-alpha_level) );
if vec_of_maxima(1) > thresh
    fpr = fpr + 1;
end

end

fpr = fpr/niters
