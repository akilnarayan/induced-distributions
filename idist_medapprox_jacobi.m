function[x0] = idist_medapprox_jacobi(alph, bet, n)

if n > 0
  x0 = (bet^2 - alph^2)./(2*n + alph + bet).^2;
else
  x0 = 2/(1 + (alph+1)/(bet+1)) - 1;
end
