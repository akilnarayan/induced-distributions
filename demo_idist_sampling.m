% Demo: Samples from an additive mixture of tensorial induced measures.
% Compares an empirical distribution from this sampling against an n-asymptotic
% conjecture for the distribution.

clear
close all

d = 2;
degree = 100;

% Number of samples to generate
N = 5e3;

% We'll sample in batches of 1e3 to avoid generating a huge set lambdas
batch_size = 1e3;
M = ceil(N/batch_size);
r = zeros([N 1]);

alph = 2;
rho = 0;

univ_inv = @(uu, nn) idistinv_freud(uu, nn, alph, rho);

for m = 1:M

  fprintf('Computing batch %d...\n', m);

  bsize = min(batch_size, N - ((m-1)*batch_size));

  i1 = (m-1)*batch_size + 1;
  i2 = i1 + bsize - 1;

  % To randomly sample indices here:
  %lambdas = sampling_total_degree_indices(bsize, d, degree);
  %x = idist_mixture_sampling(lambdas, univ_inv);

  % Or to let idist_mixture_sampling randomly subsample indices from a list:
  lambdas = total_degree_indices(d, degree);
  x = idist_mixture_sampling(bsize, lambdas, univ_inv);

  r(i1:i2) = sqrt(sum(x.^2, 2));

end

% Compare distribution functions for abs(x/sqrt(2*degree))
u = linspace(0, 1, 1e3)';

F_equil = hermite_equilibrium_distribution_conjecture(u, d);

% Create an empirical CDF for plotting
xn = sort(r)/sqrt(2*degree);
yn = [ [0; xn] ...
      [xn; max(xn(N), 1)] ];
Fyn = [ (0:N)'./N (0:N).'/N ];

yn = reshape(yn', [2*N+2 1]);
Fyn = reshape(Fyn', [2*N+2 1]);

plot(yn, Fyn, 'r', 'linewidth', 3); hold on;
plot(u, F_equil, 'k:', 'linewidth', 2);
set(xlabel('$x/\sqrt{2 n}$'), 'interpreter', 'latex');
set(ylabel('Distribution function'), 'interpreter', 'latex');
set(gca, 'fontsize', 16, 'fontweight', 'b');
axis([0 1 0 1]);
