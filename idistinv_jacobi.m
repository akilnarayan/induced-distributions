function[x] = idistinv_jacobi(u, n, alph, bet)
% [x] = idistinv_jacobi(u, n, alph, bet)
%
% Computes the inverse of the order-n induced primitive for the Jacobi
% distribution with parameters alph and bet. Uses a bisection method in
% conjunction with forward evaluation given by idist_jacobi.

assert( (all(u(:)) >= 0) && (all(u(:)) <=1 ) );
assert( (alph > -1) && (bet > -1) );
assert( all(n(:) >= 0) );

x = zeros(size(u));

supp = [-1 1];

if numel(n) == 1
  %primitive = @(x) jacobi_induced_primitive(x, n, alph, bet);
  primitive = @(xx) idist_jacobi(xx, n, alph, bet);

  % Need 2*n + K coefficients, where K is the size of the Markov-Stiltjies binning procedure
  [a,b] = jacobi_recurrence(2*n + 400, alph, bet);

  x = idist_inverse(u, n, primitive, a, b, supp);

else

  nmax = max(n(:));
  [nn, ~, bin] = histcounts(n, -0.5:(nmax+0.5));

  [a,b] = jacobi_recurrence(2*nmax + 400, alph, bet);

  for qq = 0:nmax

    flags = bin==(qq+1);

    primitive = @(xx) idist_jacobi(xx, qq, alph, bet);
    x(flags) = idist_inverse(u(flags), qq, primitive, a, b, supp);

  end

end
