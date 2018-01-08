function[data] = fidistinv_hfreud_setup(n, alph, rho, data)
% data = fidistinv_hfreud_setup(n, alph, rho, data)
%
% Setup computations for a Fast Induced Distribution Inverse routine for
% half-line Freud weights.
%
% Computes coefficients for a piecewise Chebyshev interpolant of the degree-n
% induced distribution for half-line Freud weights. These coefficients are
% appended to the data cell array, inside data{n+1}.
%
% This function also computes the "downward closed" data, that is, if
% length(data)==k, meaning data up to degree k-1 is present, then this function
% also computes and stores data for degrees k, k+1, ..., n.

ns = length(data):n;
if numel(ns) < 1
  return
end

% Number of Chebyshev coefficients per subinterval
M = 50;
tol = 1e-12;

fprintf('One-time setup computations: Computing Half-Freud (alpha=%1.3f, rho=%1.3f) induced distribution data for...\n', alph, rho);

% Construct piecewise polynomial data
for q = 1:length(ns)

  nn = ns(q);

  fprintf('n = %d...\n', nn);

  [a,b] = hfreud_recurrence(2*nn, alph, rho);
  b(1) = 1;

  xg = [0; gauss_quadrature(a, b, nn)];
  ug = idist_hfreud(xg, nn, alph, rho);

  % Also need Inf on the RHS -- just put dummy in for x for now
  ug = [ug; 1-tol];

  % Exponents for left- and right-hand sides
  exps = [rho/(rho+1), 2/3]; % the 2/3 is empirical
  [ug, exponents] = fidistinv_setup_helper1(ug, exps);
  
  % Insert some breakpoints between the penultimate and ultimate points
  M = 10; % Number of subintervals spanning last interval
  utemp = 1 - logspace(log10(1 - ug(end-1)), log10(1 - ug(end)), M+1).';
  ug = [ug(1:(end-2)); utemp];

  exponents = [exponents [zeros([1 M-1]); 2/3*ones([1 M-1])]]; % empirical

  idistinv = @(uu) idistinv_hfreud(uu, nn, alph, rho);
  M = 50; % Size of Chebyshev transform
  data{nn+1} = fidistinv_setup_helper2(ug, idistinv, exponents, M);
  
end
fprintf('Done\n');
