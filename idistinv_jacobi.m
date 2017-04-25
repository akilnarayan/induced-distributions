function[x] = jacobi_primitive_inverse(u, n, alph, bet)
% [x] = jacobi_primitive_inverse(u, n, alph, bet)
%
% Computes the inverse of the order-n induced primitive for the Jacobi
% distribution with parameters alph and bet. Uses a bisection method in
% conjunction with forward evaluation given by idist_jacobi.

assert( (all(u) >= 0) && (all(u) <=1 ) );
assert( (alph > -1) && (bet > -1) );
assert( n >= 0 );

%primitive = @(x) jacobi_induced_primitive(x, n, alph, bet);
primitive = @(x) idist_jacobi(x, n, alph, bet);
supp = [-1 1];

% Need 2*n + K coefficients, where K is the size of the Markov-Stiltjies binning procedure
[a,b] = jacobi_recurrence(2*n + 400, alph, bet);

x = idist_inverse(u, n, primitive, a, b, supp);
