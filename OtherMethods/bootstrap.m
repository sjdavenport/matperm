function bootstrap_samples = bootstrap( data, nboot, dotbootstrap, do2sample, show_loader )
% NEWFUN
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  data     a matrix of size [dim, nsubj]
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% n = 1000; data = randn(1,n); nboot = 1000;
% boot_dist = bootstrap(data, nboot, 1);
% histogram(boot_dist)
% prctile(boot_dist, 100*(1-0.025) )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if ~exist('show_loader', 'var')
    show_loader = 0;
end

if ~exist('do2sample', 'var')
    do2sample = 0;
end

if ~exist('dot', 'var')
    dotbootstrap = 0;
end
nsubj = size(data, 2);

%%  Main Function Loop
%--------------------------------------------------------------------------
bootstrap_samples = zeros(1,nboot);
data_mean = mean(data,2);

if dotbootstrap == 0
    bootstrap_samples(1) = max(sqrt(nsubj)*data_mean(:));
else
    data_std = std(data,0,2);
    bootstrap_samples(1) = max(sqrt(nsubj)*data_mean(:)./data_std(:));
end

if do2sample == 0
for I = 1:(nboot-1)
    if show_loader == 1
        loader(I-1, nboot-1, 'bootstrap progress:');
    end
    sample_mean = mean(data(:, randsample(nsubj,nsubj,1)),2);
    if dotbootstrap == 0 
        bootstrap_samples(I+1) = max(sqrt(nsubj)*(sample_mean(:) - data_mean(:)));
    else
        sample_std = std(data(:, randsample(nsubj,nsubj,1)),0, 2);
        bootstrap_samples(I+1) = max(sqrt(nsubj)*(sample_mean(:) - data_mean(:))./sample_std(:));
    end
end
elseif do2sample == 1
    for I = 1:(nboot-1)
    if show_loader == 1
        loader(I-1, nboot-1, 'bootstrap progress:');
    end
    sample_mean = mean(data(:, randsample(nsubj,nsubj,1)),2);
    if dotbootstrap == 0 
        bootstrap_samples(I+1) = max(sqrt(nsubj)*(abs(sample_mean(:) - data_mean(:))));
    else
        sample_std = std(data(:, randsample(nsubj,nsubj,1)),0, 2);
        bootstrap_samples(I+1) = max(sqrt(nsubj)*(abs(sample_mean(:) - data_mean(:)))./sample_std(:));
    end
end
end

end

