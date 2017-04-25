function[a,b] = hfreud_recurrence(N, alph, rho)
% [a,b] = hfreud_recurrence(N, alph, rho)
%
%     Returns the first N three-recurrence coefficients for the half-line Freud
%     polynomial family with parameters alph, rho.

% A limitation
if alph == 1

  N = max(N(:));
  n = (1:N).' - 1;

  a = 2*n + rho + 1;
  b = zeros(size(n));

  neq0 = (n==0);
  b(neq0) = gamma(1 + rho);
  b(~neq0) = n(~neq0).*(n(~neq0) + rho);

else 

  error('Only alpha=1 half-Freud recurrence coefficients have explicit formulas');

end
