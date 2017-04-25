function[F] = idist_freud(x, n, alph, rho)
% F = idist_freud(x, n, alph, rho)
%
%   Evaluates the integral
%
%      F = \int_{-1}^x p_n^2(x) \dx{\mu(x)},
%
%   where mu corresponds to an (alph,rho) Freud weight, scaled to be a
%   probability distribution on R, and p_n is the corresponding degree-n
%   orthonormal polynomial.

assert(rho > -1);
assert(n >= 0);

if numel(x) == 0
  F = [];
  return;
end

rflags = x > 0;
F = zeros(size(x));

F(rflags) = 1 - idist_freud(-x(rflags), n, alph, rho);

if mod(n,2) == 0
  F(~rflags) = 1/2 * idistc_hfreud(x(~rflags).^2, n/2, alph/2, (rho-1)/2);
else
  F(~rflags) = 1/2 * idistc_hfreud(x(~rflags).^2, (n-1)/2, alph/2, (rho+1)/2);
end
