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

assert( (all(u(:)) >= 0) && (all(u(:)) <=1 ) );
assert( (alph > 0) && (rho > -1) );
assert( all( n(:) >= 0 ) );

if numel(u) == 0
  x = []; 
  return
end

x = zeros(size(u));

if numel(n) == 1

  if mod(n,2) == 0
    x = sqrt( fidistinv_hfreud( abs(2*u-1), n/2, alph/2, (rho-1)/2 ) );
  else
    x = sqrt( fidistinv_hfreud( abs(2*u-1), (n-1)/2, alph/2, (rho+1)/2 ) );
  end

else

  assert(numel(n) == numel(u));

  evenflags = (mod(n,2) == 0);
  oddflags = ~evenflags;

  x(evenflags) = sqrt( fidistinv_hfreud( abs(2*u(evenflags)-1), n(evenflags)/2, alph/2, (rho-1)/2 ) );
  x(oddflags) = sqrt( fidistinv_hfreud( abs(2*u(oddflags)-1), (n(oddflags)-1)/2, alph/2, (rho+1)/2 ) );

end

flags = (u < 0.5);
x(flags) = - x(flags);
