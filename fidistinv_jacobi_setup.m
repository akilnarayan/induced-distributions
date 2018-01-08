function[data] = fidistinv_jacobi_setup(n, alph, bet, data)
% data = fidistinv_jacobi_setup(n, alph, bet, data)
%
% Setup computations for a Fast Induced Distribution Inverse routine for Jacobi
% weights.
%
% Computes coefficients for a piecewise Chebyshev interpolant of the degree-n
% induced distribution for Jacobi weights. These coefficients are appended to
% the data cell array, inside data{n+1}.
%
% This function also computes the "downward closed" data, that is, if
% length(data)==k, meaning data up to degree k-1 is present, then this function
% also computes and stores data for degrees k, k+1, ..., n.

ns = length(data):n;
if numel(ns) < 1
  return
end

fprintf('One-time setup computations: Computing Jacobi (alpha=%1.3f, beta=%1.3f) induced distribution data for...\n', alph, bet);

% Construct piecewise polynomial data
for q = 1:length(ns)

  nn = ns(q);

  fprintf('n = %d...\n', nn);

  [a,b] = jacobi_recurrence(nn+1, alph, bet);
  xg = gauss_quadrature(a,b,nn);
  ug = idist_jacobi(xg, nn, alph, bet, 50); % Make it very accurate

  ug = [0;  ug;  1];
  exps = [bet/(bet+1) alph/(alph+1)];
  [ug, exponents] = fidistinv_setup_helper1(ug, exps);

  idistinv = @(uu) idistinv_jacobi(uu, nn, alph, bet);
  M = 50; % Size of Chebyshev transform
  data{nn+1} = fidistinv_setup_helper2(ug, idistinv, exponents, M);

end
fprintf('Done\n');
