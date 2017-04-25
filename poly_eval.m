function[p] = poly_eval(a, b, x, N)
% poly_eval -- Evaluates orthogonal polynomials
%
% p = poly_eval(a, b, x, N)
%   Uses the recurrence coefficients a and b to evaluate p_n(x), 
%   where p_n(x) is the degree-n orthonormal polynomial associated with the
%   recurrence coefficients a, b (with positive leading coefficient).
%
%   The p_n satisfy the recurrences
%
%     sqrt(b_{n+1}) p_{n+1} = (x - a_n) p_n - sqrt(b_n) p_{n-1}
%
%   With the input arrays a and b, we have a_n = a(n+1), b_n = b(n+1). Note
%   that b_0 = b(1) is only used to initialize p_0.
%
%   The output matrix p has size numel(x) x N+1, and hence the first N+1 (up to
%   degree N) polynomials are evaluated.
%
%   Inputs:
%       x : array of doubles
%       N : positive integer, N - 1 < length(a) == length(b)
%       a : array of recurrence coefficients
%       b : array of reucrrence coefficients

nx = numel(x);

assert(N >= 0);

assert(N <= length(a));
assert(N <= length(b));

p = zeros([nx N+1]);

% Flatten x
xf = x(:);

p(:,1) = 1/sqrt(b(1)) * ones([nx 1]);
if N > 0
  p(:,2) = 1/sqrt(b(2)) * (xf - a(1)).*p(:,1);
end

for q = 2:N
  % Derived from three-term recurrence
  p(:,q+1) = (xf - a(q)).*p(:,q) - sqrt(b(q)).*p(:,q-1);
  p(:,q+1) = 1/sqrt(b(q+1)) * p(:,q+1);
end
