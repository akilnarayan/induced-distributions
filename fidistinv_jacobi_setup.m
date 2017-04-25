function[data] = fidistinv_jacobi(n, alph, bet, data)
% data = fidistinv_jacobi(n, alph, bet, data)
%
% Setup computations for a Fast Induced Distribution Inverse routine for Jacobi
% weights.
%
% Computes coefficients for a piecewise cubic interpolant of the degree-n
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

% Number of points per sub interval in x
Q = 20;
xtemplate = linspace(0, 1, Q); xtemplate(1) = [];
xtemplate = (sort(cos(pi*xtemplate)) + 1)/2;

% Enrich with R points equidistantly in both x and u 
R = 100;

fprintf('One-time setup computations: Computing induced distribution data for...\n');

% Construct piecewise polynomial data
for q = 1:length(ns)

  nn = ns(q);

  fprintf('n = %d...\n', nn);
  
  % The data points we'll choose are as follows:
  % R equally spaced data on u and x space
  % Q data points on each interval bounded by zeros of p_nn (Q*(nn+1) + 1 points)

  [a,b] = jacobi_recurrence(nn+1, alph, bet);
  x = gauss_quadrature(a, b, nn);

  x = [-1; x; 1];

  % Create Q-1 points per interval
  xn = diff(x) * xtemplate + repmat(x(1:end-1), [1 Q-1]);
  xn = [reshape(xn.', [numel(xn) 1]); 1];
  X = numel(xn);

  un = zeros(size(xn));
  un = idist_jacobi(xn, nn, alph, bet);

  % Now R data points equally-spaced in u space
  temp = linspace(0, 1, R + 2).';
  temp([1 end]) = [];

  un = [un; temp];
  xn = [xn; idistinv_jacobi(temp, nn, alph, bet)];

  % Now R data points equally-spaced in x space
  temp = linspace(-1, 1, R+2).';
  temp([1 end]) = [];

  X = numel(xn);
  xn = [xn; temp];
  un = [un; zeros(size(temp))];
  un((X+1):numel(xn)) = idist_jacobi(xn((X+1):numel(xn)), nn, alph, bet);

  % Sort
  [xn, inds] = sort(xn);
  un = un(inds);

  assert( all(diff(un) >= 0), 'Induced function values are not monotonic' );

  data{nn+1} = pchip(xn, un);
  data{nn+1}.u = un;

end
fprintf('Done\n');
