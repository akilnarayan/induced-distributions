function[x] = fidistinv_jacobi(u, n, alph, bet);
% [x] = fidistinv_jacobi(u, n, alph, bet)
%
% A Fast Induced Distribution Inverse routine for Jacobi weights.
%
% Computes the inverse of the order-n induced primitive for the Jacobi
% distribution with parameters alph and bet. A monotone piecewise cubic
% interpolant of the distribution function is computed and saved for future
% runs.

data = load_fjacobi(n, alph, bet);

if length(data) < max(n(:))+1

  data = fidistinv_jacobi_setup(max(n(:)), alph, bet, data);
  save_fjacobi(data, alph, bet);

end

x = fidistinv_driver(u, n, data);
