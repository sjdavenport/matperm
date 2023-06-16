spfn = @(n) randn(n, 1)'; nsubj = 1000;
alpha = 0.05; niters = 1000;

gamma = [0.4:0.1:0.9];
fpr_blb = zeros(1, length(gamma));
for I = 1:length(gamma)
    I
    fpr_blb(I) = rc_blb(  spfn, nsubj, floor(nsubj^(gamma(I))), 5, 30, niters, alpha );
end

%%
% gamma = [0.2:0.1:1];
niters = 1000;
% load('./onesamplepower.mat')
plot(gamma, fpr_blb, 'color', 'blue');
std_error = bernstd( 0.05, niters, 0.95 );
hold on
plot(gamma, std_error(1)*ones(1,length(gamma))', '--', 'color', 'blue');
plot(gamma, std_error(2)*ones(1,length(gamma))', '--', 'color', 'blue');

%%
spfn = @(n) randn(n, 1)' + 0.01; nsubj = 1000;
alpha = 0.05; niters = 1000;

fpr_blb = rc_blb(  spfn, nsubj, 100, 5, 30, niters, alpha )

%%
spfn = @(n) randn(n, 1)' + 0.01; nsubj = 1000;
alpha = 0.05; niters = 1000;

fpr_blb = rc_blb(  spfn, nsubj, 100, 5, 30, niters, alpha )
