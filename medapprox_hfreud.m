function[x0] = medapprox_hfreud(alph, rho, n)
% x0 = medapprox_hfreud(alph, rho, n)
%
% Returns a guess for the median of the order-n half-Freud induced distribution
% with parameters alph, rho.

x0 = 1/2 * ( maxapprox_hfreud(alph, rho, n) + minapprox_hfreud(alph, rho, n) );
