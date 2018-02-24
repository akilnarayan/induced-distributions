function[x] = idist_inverse(u, n, primitive, a, b, supp)
% [x] = idist_inverse(u, n, primitive, a, b, supp)
%
% Uses bisection to compute the (approximate) inverse of the order-n induced
% primitive function F_n. 
%
% The ouptut x = F_n^{-1}(u). The input function primitive should be a function
% handle accepting a single input and outputs the primitive F_n evaluated at
% the input. 
%
% The last inputs a, b, and supp are three-term recurrence coefficients (a and b) for the
% original measure, and its support (supp). These are required to formulate a
% good initial guess via the Markov-Stiltjies inequalities.

persistent options
if isempty(options)
  %% fzero
  options = optimset(@fzero);
  options.Display = 'off';
  options.TolFun = 10*eps;
  options.TolX = 10*eps;
end

if numel(n) == 1
  intervals = markov_stiltjies_initial_guess(u, n, a, b, supp);
else
  intervals = zeros([numel(n) 2]);
  nmax = max(n(:));
  [nn, ~, bin] = histcounts(n, -0.5:(0.5+nmax));
  for qq = 0:max(n(:))
    flags = bin==(qq+1);
    intervals(flags,:) = markov_stiltjies_initial_guess(u(flags), qq, a, b, supp);
  end
end

x = zeros(size(u));

for q = 1:numel(u)
  fun = @(xx) primitive(xx) - u(q);
  x(q) = fzero(fun, intervals(q,:), options);
end
