function[intervals] = markov_stiltjies_initial_guess(u, n, a, b, supp)
% intervals = markov_stiltjies_initial_guess(u, n, a, b, supp)
%
% Uses the Markov-Stiltjies inequalities to provide a bounding interval for x
% the solution to 
%
%   F_n(x) = u,
%
% where n is the the order-n induced distribution function associated to the
% measure with three-term recurrrence coefficients a, b, having support on the
% real-line interval defined by the length-2 vector supp.
%
% If u is a length-M vector, the output intervals is an (M x 2) matrix, with
% row m the bounding interval for u = u(m).

assert(numel(a) == numel(b));
assert(numel(a) > 2*max(n(:)));

% Compute quadratic modifications modifications.
[x,w] = gauss_quadrature(a, b, n);
b(1) = 1;
for k = 1:n
  [a,b] = quadratic_modification_C(a, b, x(k));
  b(1) = 1;
end

%% Markov-Stiltjies inequalities
% Use all the remaining coefficients for the Markov-Stiltjies inequalities
N = length(a);
[y,w] = gauss_quadrature(a, b, N);
if supp(2) > y(end)
  X = [supp(1); y; supp(2)];
  W = [0; cumsum(w)]; 
else
  X = [supp(1); y; y(end)];
  W = [0; cumsum(w)]; 
end
W = W/W(end);

W(W > 1) = 1;% Just in case for machine eps issues
W(end) = 1; 

[~,j] = histc(u, W);
j = j(:);
jleft = j;
jright = jleft + 2;

% Fix endpoints
flags = (jleft == (N+1));
jleft(flags) = N+2;
jright(flags) = N+2;

intervals = [X(jleft) X(jright)];
