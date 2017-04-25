% Demo: computes error in M-point quadrature rule for Jacobi weights.

clear
close all

% The "exact" answer
M_exact = 1e3;

a = exp(1);
b = -1/3;
n = 2;

xs = linspace(-1, 1, 201).';

Ms = (1:10).';

exact = zeros([1 numel(xs)]);
F = zeros([numel(xs) numel(Ms)]);

errors = zeros([numel(xs) numel(Ms)]);

exact = idist_jacobi(xs, n, a, b, M_exact);
for q_m = 1:length(Ms)
  F(:,q_m) = idist_jacobi(xs, n, a, b, Ms(q_m));
  fprintf('M = %d\n', Ms(q_m));
end

errors = abs(F - repmat(exact, [1 length(Ms)]));

[xx, MM] = ndgrid(xs, Ms);

logerrs = log10(errors');
logerrs(isinf(logerrs)) = -16;
h = pcolor(xx', MM', logerrs);
set(h, 'EdgeColor','none', 'xdata', xs.');
shading interp;
set(colorbar, 'fontsize', 16, 'fontweight', 'b');
caxis([-16 0]);
set(gca, 'fontsize', 16, 'fontweight', 'b');
h = ylabel('$M$');
set(ylabel('$\boldmath{M}$'), 'interpreter', 'latex', 'fontsize', 16);
set(xlabel('$\boldmath{x}$'), 'interpreter', 'latex', 'fontsize', 16);
colormap hot

hold on;
med = medapprox_jacobi(a, b, n);
plot([med med], Ms([1 end]), 'k--', 'linewidth', 2);
