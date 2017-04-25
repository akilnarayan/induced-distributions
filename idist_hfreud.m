function[F] = idist_hfreud(x, n, alph, rho, M)
% idist_hfreud -- Evaluation of induced distribution
%
% F = idist_hfreud(x, n, alph, rho, M)
%
%   Evaluates the integral
%
%      F = \int_{0}^x p_n^2(x) \dx{\mu(x)},
%
%   where mu corresponds to a (alph,rho) generalized Freud weight, scaled to be a
%   probability distribution on R, and p_n is the corresponding degree-n
%   orthonormal polynomial.
%
%   This function evaluates this via a transformation, measure modification,
%   and Gauss quadrature, the ending Gauss quadrature has M points.

assert(rho > -1);
assert( all(x >= 0) );
assert(n >= 0);

if numel(x) == 0
  F = [];
  return;
else
  F = zeros(size(x));
end

if nargin < 5
  if alph ~= 1
    M = n + 10; % Complementary integral is much less accurate for alph ~= 1
  else
    M = 25;
  end
end

if alph ~= 1

  x0 = medapprox_hfreud(alph, rho, n);
  rflags = x > x0;
  F(rflags) = 1 - idistc_hfreud(x(rflags), n, alph, rho, M);

else

  x0 = 50;
  rflags = x > x0;

  F(rflags) = 1 - idistc_hfreud(x(rflags), n, alph, rho, M);

end

if n > 0

  [a,b] = hfreud_recurrence(n+1, alph, rho);
  b(1) = 1;

  % Zeros of p_n
  xn = gauss_quadrature(a, b, n);

  % This is the (inverse) n'th root of the leading coefficient square of p_n
  % We'll use it for scaling later
  logfactor0 = -sum(log(b));
else
  logfactor0 = 0;
  xn = [];
end

% Jacobi recurrence coefficients
[aJ, bJ] = jacobi_recurrence(2*n+M+1, 0, rho);
bJ(1) = 1;

for xq = 1:numel(x)

  if x(xq) == 0
    F(xq) = 0;
    continue
  end

  if rflags(xq)
    continue
  end

  % Transformed zeros
  un = 2*xn/x(xq) - 1; 

  a = aJ; 
  b = bJ; 
  logfactor = logfactor0;

  % Successive quadratic measure modifications
  for j = 1:n
    [a,b] = quadratic_modification_C(a, b, un(j));

    logfactor = logfactor + log(b(1));
    b(1) = 1;
  end

  % M-point Gauss quadrature for evaluation of auxilliary integral I
  [u,w] = gauss_quadrature(a, b, M);
  I = w.'*exp(- (x(xq)/2)^alph * (u+1).^alph);

  logfactor = logfactor + (2*n+rho+1)*log(x(xq)/2) + log(alph) + (rho+1)*log(2) - gammaln((rho+1)/alph) - log(rho+1);

  F(xq) = exp(logfactor + log(I));

end
