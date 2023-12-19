function [ threshold, vec_of_max_medians, permuted_mean_store ] = ... 
        perm_mean( data, alpha, nperm, show_loader, store_perms )
% NEWFUN
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  data: 
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% dim = [50,50]; nsubj = 50; FWHM = 0;
% Sig = 0.5*peakgen(1, 10, 8, dim);
% data = wfield(dim, nsubj).field + Sig;
% data_vec = vec_data(data, ones(dim));
% threshold = perm_median(data_vec)
% data_median = median(data, 3);
% sum(data_median(:) > threshold)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
s_data = size(data);
nsubj = s_data(end);

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'nperm', 'var' )
   % Default value
   nperm = 1000;
end

if ~exist( 'show_loader', 'var' )
   % Default value
   show_loader = 1;
end

if ~exist( 'alpha', 'var' )
   % Default value
   alpha = 0.05;
end

if ~exist('store_perms', 'var')
    store_perms = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
mtstat = mean(data,2);
max_mean = max(mtstat);

vec_of_max_medians = zeros(1,nperm);
vec_of_max_medians(1) = max_mean;

% Compute bernoulli random variables for the sign flipping
random_berns = 2*(binornd(1,0.5, nsubj, nperm )-1/2);

if store_perms == 1
    permuted_mean_store = zeros([s_data(1), nperm]);
    permuted_mean_store(:,1) = max_mean;
else
    permuted_mean_store = NaN;
end

for I = 2:nperm
    if show_loader == 1
        loader(I-1, nperm-1, 'perm progress:');
    end
    
    random_berns_for_iter = random_berns(:, I);
    random_sample_negative = find(random_berns_for_iter < 0);

    data_perm = data;
    data_perm(:,random_sample_negative) = -data(:,random_sample_negative);
    
    mstat_perm = mean(data_perm,2);
    max_mean_perm = max(mstat_perm);
    
    if store_perms == 1
        permuted_mean_store(:,I) = max_mean_perm;
    end
    
    vec_of_max_medians(I) = max_mean_perm;
end

threshold = prctile(vec_of_max_medians, 100*(1-alpha));

end

