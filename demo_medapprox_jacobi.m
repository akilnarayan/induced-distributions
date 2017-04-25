% Demo: compares actual induced distribution medians versus the
% potential-theroetic guesses for Jacobi weights.

clear
close all

ns = 1:300;

alph1 = 385;
bet1 = pi;

alph2 = -1/pi;
bet2 = pi*100;

alph3 = 0.5;
bet3 = 2;

potential_theory_guesses1 = medapprox_jacobi(alph1, bet1, ns);
potential_theory_guesses2 = medapprox_jacobi(alph2, bet2, ns);
potential_theory_guesses3 = medapprox_jacobi(alph3, bet3, ns);

% Now compute medians
medians1 = zeros(size(ns));
medians2 = zeros(size(ns));
medians3 = zeros(size(ns));
for q = 1:length(ns);

  n = ns(q);

  medians1(q) = idistinv_jacobi(0.5, n, alph1, bet1);
  medians2(q) = idistinv_jacobi(0.5, n, alph2, bet2);
  medians3(q) = idistinv_jacobi(0.5, n, alph3, bet3);

  if mod(q,10)==0
    fprintf('n = %d\n', n);
  end

end

plot(potential_theory_guesses1, ns, 'r', medians1, ns, 'r.');
hold on;
plot(potential_theory_guesses2, ns, 'b', medians2, ns, 'b.');
plot(potential_theory_guesses3, ns, 'k', medians3, ns, 'k.');
axis([-1 1 1 300]);
set(xlabel('$x$'), 'interpreter', 'latex', 'fontsize', 16, 'fontweight', 'b');
set(ylabel('$n$'), 'interpreter', 'latex', 'fontsize', 16, 'fontweight', 'b');
set(gca, 'fontsize', 16, 'fontweight', 'b');
