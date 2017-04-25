function[x0] = medapprox_jacobi(alph, bet, n)
% x0 = medapprox_jacobi(alph, bet, n)
%
% Returns a guess for the median of the order-n Jacobi induced distribution
% with parameters alph, bet.

if n > 0
  x0 = (bet^2 - alph^2)./(2*n + alph + bet).^2;
else
  x0 = 2/(1 + (alph+1)/(bet+1)) - 1;
end
