clear
close all

% The "exact" answer
M_exact = 1000;

alph = 1;
rho = sqrt(1001);
n = 595;

xs = linspace(0, 100, 51)';

Ms = (10:25).';

exact = zeros([1 numel(xs)]);
F = zeros([numel(xs) numel(Ms)]);

errors = zeros([numel(xs) numel(Ms)]);

fprintf('Computing ''exact'' answer with M = %d...\n', M_exact);
exact = idist_hfreud(xs, n, alph, rho, M_exact);

for q_m = 1:length(Ms)
  fprintf('M = %d...\n', Ms(q_m));
  F(:,q_m) = idist_hfreud(xs, n, alph, rho, Ms(q_m));
  errors(:,q_m) = abs(F(:,q_m) - exact);
end

[xx, MM] = ndgrid(xs, Ms);

h = pcolor(xx', MM', log10(errors'));
set(h, 'EdgeColor','none', 'xdata', xs.');
shading interp;
set(colorbar, 'fontsize', 16, 'fontweight', 'b');
set(gca, 'fontsize', 16, 'fontweight', 'b');
h = ylabel('$M$');
set(ylabel('$\boldmath{M}$'), 'interpreter', 'latex', 'fontsize', 16);
set(xlabel('$\boldmath{x}$'), 'interpreter', 'latex', 'fontsize', 16);
caxis([-16 0]);

hold on;
med = medapprox_hfreud(alph, rho, n);
plot([med med], Ms([1 end]), 'k');
colormap('hot');
