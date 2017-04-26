function[x] = fidistinv_freud(u, n, alph, rho)
% [x] = fidistinv_freud(u, n, alph, rho)
%
% A Fast Induced Distribution Inverse routine for Freud weights.
%
% Computes the inverse of the order-n induced primitive for the Freud
% distribution with parameters alph and rho. Uses the methodology in
% fidistinv_hfreud.m.
%
% Supports vectorization in u.

assert( (all(u) >= 0) && (all(u) <=1 ) );
assert( (alph > 0) && (rho > -1) );
assert( n >= 0 );

if numel(u) == 0
  x = []; 
  return
end

x = zeros(size(u));

if mod(n,2) == 0
  x = sqrt( fidistinv_hfreud( abs(2*u-1), n/2, alph/2, (rho-1)/2 ) );
else
  x = sqrt( fidistinv_hfreud( abs(2*u-1), (n-1)/2, alph/2, (rho+1)/2 ) );
end
flags = (u < 0.5);
x(flags) = - x(flags);
