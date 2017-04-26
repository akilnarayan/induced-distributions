% Demo: illustrates sampling via fast and standard methods

clear
close all

alph = 3;
bet = -0.3;
n = 10;

N = 1e3;

u = rand([N 1]);

% The standard way
tic;
x = idistinv_jacobi(u, n, alph, bet);
time = toc/N;

% The fast way
tic
xf = fidistinv_jacobi(u, n, alph, bet);
ftime = toc/N;

% Multivariate fast sampling of additive mixture of induced distributions
d = 2;
lambdas = total_degree_indices(d, n);

% The original way
univ_inv = @(uu,nn) idistinv_jacobi(uu, nn, alph, bet);
tic;
x2 = idist_sampling(N, lambdas, univ_inv);
time2 = toc/N;

% The fast way
univ_inv = @(uu,nn) fidistinv_jacobi(uu, nn, alph, bet);
tic;
xf2 = idist_sampling(N, lambdas, univ_inv);
ftime2 = toc/N;
