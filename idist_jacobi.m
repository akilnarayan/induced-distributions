function[F] = idist_jacobi(x, n, alph, bet, M)
% idist_jacobi -- Evaluation of induced distribution
%
% F = idist_jacobi(x, n, alph, bet, {M = 10})
%
%   Evaluates the integral
%
%      F = \int_{-1}^x p_n^2(x) \dx{\mu(x)},
%
%   where mu is the (a,b) Jacobi polynomial measure, scaled to be a probability
%   distribution on [-1,1], and p_n is the corresponding degree-n orthonormal
%   polynomial.
%
%   This function evaluates this via a transformation, measure modification,
%   and Gauss quadrature, the ending Gauss quadrature has M points.

assert((alph > -1) && (bet > -1));
assert( all(abs(x) <= 1) );
assert( (n >= 0) && (numel(n) == 1) );

if numel(x) == 0
  F = [];
  return
end

if nargin < 5
  M = 10;
end

A = floor(abs(alph)); % is an integer
Aa = alph - A;

F = zeros(size(x));

mrs_centroid = medapprox_jacobi(alph, bet, n);
xreflect = x > mrs_centroid;

F(xreflect) = 1 - idist_jacobi(-x(xreflect), n, bet, alph, M);

[a,b] = jacobi_recurrence(n+1, alph, bet);
b(1) = 1; % To make it a probability measure

if n > 0
  % Zeros of p_n
  xn = gauss_quadrature(a, b, n);
end

% This is the (inverse) n'th root of the leading coefficient square of p_n
% We'll use it for scaling later
kn_factor = exp(-1/n*sum(log(b)));

for xq = 1:numel(x)

  if x(xq) == -1
    F(xq) = 0;
    continue
  end

  if xreflect(xq)
    continue
  end

  % Recurrence coefficients for quadrature rule
  [a,b] = jacobi_recurrence(2*n+A+M+1, 0, bet);
  b(1) = 1; % To make it a probability measure
  a = a.'; b = b.';

  if n > 0
    % Transformed
    un = (2/(x(xq)+1)) * (xn + 1) - 1;
  end

  logfactor = 0; % Keep this so that bet(1) always equals what it did before

  % Successive quadratic measure modifications
  for j = 1:n
    [a,b] = quadratic_modification_C(a, b, un(j));

    logfactor = logfactor + log(b(1)*((x(xq)+1)/2)^2 * kn_factor);
    b(1) = 1;
  end

  % Linear modification by factors (2 - 1/2*(u+1)*(x+1)), having root u = (3-x)/(1+x)
  root = (3-x(xq))/(1+x(xq));
  for aq = 1:A
    [a, b] = linear_modification(a, b, root);

    logfactor = logfactor + log(b(1) * 1/2 * (x(xq)+1));
    b(1) = 1;
  end

  % M-point Gauss quadrature for evaluation of auxilliary integral I
  [u,w] = gauss_quadrature(a, b, M);
  I = w.'*(2 - 1/2 * (u+1) * (x(xq)+1) ).^Aa;

  F(xq) = exp(logfactor - alph*log(2) - betaln(bet+1, alph+1) - log(bet+1) + (bet+1)*log((x(xq)+1)/2)) * I;

end
