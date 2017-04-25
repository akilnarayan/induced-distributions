function[C] = C_eval(a, b, x, N)
% C_eval -- Evaluates Christoffel-normalized orthogonal polynomials
%
% C = C_eval(a, b, x, N)
%   Uses the recurrence coefficients a and b to evaluate C_n(x), 
%   defined as
%
%    C_n(x) = p_n(x) / sqrt(sum_{j=0}^{n-1} p_j^2(x)),
%
%   where p_n(x) is the degree-n orthonormal polynomial associated with the
%   recurrence coefficients a, b (with positive leading coefficient). We define
%   C_{-1} = 0, and C_0 = p_0.
%
%   The p_n satisfy the recurrences
%
%     sqrt(b_{n+1}) p_{n+1} = (x - a_n) p_n - sqrt(b_n) p_{n-1}
%
%   With the input arrays a and b, we have a_n = a(n+1), b_n = b(n+1). Note
%   that b_0 = b(1) is only used to initialize p_0.
%
%   The output matrix C has size numel(x) x N, and hence the first N (up to
%   degree N-1) polynomials are evaluated.
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

C = zeros([nx N+1]);

% Flatten x
xf = x(:);

% To initialize C, we need p_0 and p_1
C(:,1) = 1./sqrt(b(1)); % C0 = p0
if N > 0
  C(:,2) = 1/sqrt(b(2)) * (xf - a(1));
  % C1 = p1/p0
end

if N > 1
  C(:,3) = 1./sqrt(1 + C(:,2).^2) .* ( (xf - a(2)).*C(:,2) - sqrt(b(2)) );
  C(:,3) = C(:,3)./sqrt(b(3));
end

for n = 3:N
  C(:,n+1) = 1./sqrt(1 + C(:,n).^2) .* ( (xf - a(n)).*C(:,n) - sqrt(b(n))*C(:,n-1)./sqrt(1 + C(:,n-1).^2));
  C(:,n+1) = C(:,n+1)/sqrt(b(n+1));
end
