function[x] = hfreud_tolerance(n, alph, rho, tol)
% [x] = hfreud_tolerance(n, alph, rho, tol)
%
% Computes a value x such that
%
%     F_n(x) >= 1 - tol,
%
% where F_n is the order-n induced distribution for the half-Freud measure with
% parameters (alph, rho).

assert( (0 < tol) && (tol < 1) );

x = maxapprox_hfreud(alph, rho, n);
F = idistc_hfreud(x, n, alph, rho);

while F > tol

  x = x + 1;
  F = idistc_hfreud(x, n, alph, rho);

end
