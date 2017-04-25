% Demo: compares potential-theoretic guesses for the essential support interval
% of the order-n induced distribution versus actual support intervals.

clear
close all

alph = 1;
rho = 20*pi;

ns = 0:50;

% Compute bounding intervals:
intervals = zeros([numel(ns) 2]);
Fbounds = zeros([numel(ns) 2]);

for q = 1:numel(ns)
  n = ns(q);

  intervals(q,1) = minapprox_hfreud(alph, rho, n);
  intervals(q,2) = maxapprox_hfreud(alph, rho, n);

  Fbounds(q,1) = idistinv_hfreud(0.01, n, alph, rho);
  Fbounds(q,2) = idistinv_hfreud(0.99, n, alph, rho);

  if mod(q, 10) == 0
    fprintf('n = %d\n', n);
  end
end

figure;
plot(intervals(:,1), ns, 'r', 'linewidth', 2); hold on 
plot(intervals(:,2), ns, 'r', 'linewidth', 2);
for q = 1:numel(ns)
  n = ns(q);
  plot(Fbounds(q,:), [n n], 'k', 'linewidth', 2);
end
set(xlabel('$x$'), 'interpreter', 'latex');
set(ylabel('$n$'), 'interpreter', 'latex');
set(gca, 'fontsize', 16, 'fontweight', 'b');

xmax = max(max(intervals(:,2)), max(Fbounds(:,2)));
axis([0, 1.05*xmax, min(ns) - 1, max(ns)+1]);
