% Demo: illustrates sampling via fast and standard methods

clear
close all

alph = 1;
rho = pi;
n = 10;

N = 1e3;

u = rand([N 1]);

% The standard way
tic;
x = idistinv_hfreud(u, n, alph, rho);
time = toc/N;

% The fast way
tic
xf = fidistinv_hfreud(u, n, alph, rho);
ftime = toc/N;

% Multivariate fast sampling of additive mixture of induced distributions
d = 2;
lambdas = total_degree_indices(d, n);

% The original way
univ_inv = @(uu,nn) idistinv_hfreud(uu, nn, alph, rho);
tic;
x2 = idist_mixture_sampling(N, lambdas, univ_inv);
time2 = toc/N;

% The fast way
univ_inv = @(uu,nn) fidistinv_hfreud(uu, nn, alph, rho);
tic;
xf2 = idist_mixture_sampling(N, lambdas, univ_inv);
ftime2 = toc/N;
