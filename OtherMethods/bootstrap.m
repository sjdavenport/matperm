function bootstrap_samples = bootstrap( data, nboot )
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
% n = 10000; data = randn(1,n); nboot = 1000;
% boot_dist = bootstrap(data, nboot);
% histogram(boot_dist)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

nsubj = size(data, 2);

%%  Main Function Loop
%--------------------------------------------------------------------------
bootstrap_samples = zeros(1,nboot);
data_mean = mean(data,2);
bootstrap_samples(1) = max(sqrt(nsubj)*data_mean(:));

for I = 1:(nboot-1)
    sample_mean = mean(data(:, randsample(nsubj,nsubj,1)),2);
    bootstrap_samples(I+1) = max(sqrt(nsubj)*(sample_mean(:) - data_mean(:)));
end

end

