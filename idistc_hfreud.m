function[F] = idistc_hfreud(x, n, alph, rho, M)
% idistc_hfreud -- Evaluation of complementary induced distribution
%
% F = idistc_hfreud(x, n, alph, rho)
%
%   Evaluates the integral
%
%      F = \int_{x}^\infty p_n^2(x) \dx{\mu(x)},
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
  lflags = x <= x0;
  F(lflags) = 1 - idist_hfreud(x(lflags), n, alph, rho, M);

else

  x0 = 50;
  lflags = x <= x0;

  F(lflags) = 1 - idist_hfreud(x(lflags), n, alph, rho, M);

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

% Compute measure modification associated to rho
R = floor(abs(rho));

% Compute recurrence coefficients
[aH,bH] = hfreud_recurrence(2*n+M+R+1, alph, 0);
bH(1) = 1; % To make it a probability measure

for xq = 1:numel(x)

  if x(xq) == Inf
    F(xq) = 0;
    continue
  end

  if lflags(xq)
    continue
  end

  % Transformed
  un = xn - x(xq);

  a = aH;
  b = bH;
  logfactor = logfactor0;

  % Successive quadratic measure modifications
  for j = 1:n
    [a,b] = quadratic_modification_C(a, b, un(j));

    logfactor = logfactor + log(b(1));
    b(1) = 1;
  end

  root = -x(xq);
  for rq = 1:R
    [a,b] = linear_modification(a, b, root);

    logfactor = logfactor + log(b(1));
    b(1) = 1;
  end

  % M-point Gauss quadrature for evaluation of auxilliary integral I
  [u,w] = gauss_quadrature(a, b, M);
  I = w.'*((u + x(xq)).^(rho-R).*exp(u.^alph + x(xq)^alph - (u+x(xq)).^alph));

  % Instead of this, we build this dependence into kn_factor
  F(xq) = I* exp(-x(xq)^alph + logfactor + gammaln(1/alph) - gammaln((rho+1)/alph));

end
