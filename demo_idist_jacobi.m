clear
close all

% Plots showing the induced distribution and density for a certain (alph,bet,n)
% Jacobi case.

alph = -0.8;
bet = sqrt(101);
n = 8;

M = 1001;

x = cos(linspace(0, pi, M+2).'); x([1 end]) = [];

F = idist_jacobi(x, n, alph, bet);

[a,b] = jacobi_recurrence(n+1, alph, bet); b(1) = 1;
polys = poly_eval(a, b, x, n);
f = polys(:,end).^2.*(1-x).^alph.*(1+x).^bet;
f = f/ (2^(alph+bet+1) * beta(bet+1, alph+1));

figure;
plot(x, f);
set(xlabel('$x$'), 'interpreter', 'latex');
set(ylabel('$p_n^2(x) \mathrm{d}\mu(x)$'), 'interpreter', 'latex');
set(gca, 'fontsize', 16, 'fontweight', 'b');
axis([-1 1 0 4]);

figure;
plot(x, F);
set(xlabel('$x$'), 'interpreter', 'latex');
set(ylabel('$F_n(x)$'), 'interpreter', 'latex');
set(gca, 'fontsize', 16, 'fontweight', 'b');
