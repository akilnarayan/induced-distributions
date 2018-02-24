function[x] = idistinv_hfreud(u, n, alph, rho)
% [x] = idistinv_hfreud(u, n, alph, rho)
%
% Computes the inverse of the order-n induced primitive for the Half-line Freud
% distribution with parameters alph and rho. Uses a bisection method in
% conjunction with forward evaluation given by idist_hfreud.

assert( (all(u(:)) >= 0) && (all(u(:)) <=1 ) );
assert( (alph > 0) && (rho > -1) );
assert( all(n(:) >= 0) );

if numel(n) == 1
  rhs = 1.2*maxapprox_hfreud(alph, rho, n);

  U = max(u(:));
  if U == 1
    rhs = hfreud_tolerance(n, alph, rho, eps/10);
  else
    rhs = hfreud_tolerance(n, alph, rho, 1 - U);
  end

  supp = [0 rhs];

  % Need 2*n + K coefficients, where K is the size of the Markov-Stiltjies binning procedure
  [a,b] = hfreud_recurrence(2*n + max(100,n), alph, rho);

  primitive = @(x) idist_hfreud(x, n, alph, rho);
  x = idist_inverse(u, n, primitive, a, b, supp);

else

  nmax = max(n(:));
  [nn, ~, bin] = histcounts(n, -0.5:(nmax+0.5));
  [a,b] = hfreud_recurrence(2*n + max(100,nmax), alph, rho);

  for qq = 0:nmax

    rhs = 1.2*maxapprox_hfreud(alph, rho, qq);

    U = max(u(:));
    if U == 1
      rhs = hfreud_tolerance(qq, alph, rho, eps/10);
    else
      rhs = hfreud_tolerance(qq, alph, rho, 1 - U);
    end

    supp = [0 rhs];

    flags = bin==(qq+1);

    primitive = @(x) idist_hfreud(x, qq, alph, rho);

    x(flags) = idist_inverse(u(flags), qq, primitive, a, b, supp);

  end

end

