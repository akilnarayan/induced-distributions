function[a] = maxapprox_hfreud(alph, rho, n)
% a = maxapprox_hfreud(alph, rho, n)
%
% Returns a guess for the right-hand side of the support interval for
%
%     x^rho p_{n}^2(x) exp(-x^alph)
%
% where p_n is the family of orthonormal polynomials for the half-Freud weight
% weight parameters alph, rho.

if n > 0

  a = rho + 2*n + 2*sqrt(rho*n + n.^2);
  a = a.^(1/alph);
  a = a * exp(1/alph * (log(sqrt(pi)) - log(2) + gammaln(alph) - gammaln(alph+1/2)));

else

  a = gammaincinv( (1 - 1e-3) , (rho+1)/alph);
  a = a.^(1/alph);

end
