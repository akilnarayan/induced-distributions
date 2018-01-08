function[x] = fidistinv_hfreud(u, n, alph, rho)
% [x] = fidistinv_hfreud(u, n, alph, rho)
%
% A Fast Induced Distribution Inverse routine for Half-line Freud weights.
%
% Computes the inverse of the order-n induced primitive for the half-line
% Freud distribution with parameters alph and rho.

data = load_fhfreud(n, alph, rho);

if length(data) < max(n(:))+1

  data = fidistinv_hfreud_setup(max(n(:)), alph, rho, data);
  save_fhfreud(data, alph, rho);

end

x = fidistinv_driver(u, n, data);
