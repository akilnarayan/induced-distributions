function[x] = fidistinv_jacobi(u, n, alph, bet)
% [x] = fidistinv_jacobi(u, n, alph, bet)
%
% A Fast Induced Distribution Inverse routine for Jacobi weights.
%
% Computes the inverse of the order-n induced primitive for the Jacobi
% distribution with parameters alph and bet. A monotone piecewise cubic
% interpolant of the distribution function is computed and saved for future
% runs.

persistent options
if isempty(options)
  options = optimset(@fzero);
  options.Display = 'off';
  options.TolFun = 10*eps;
  options.TolX = 10*eps;
end

data = load_fjacobi(n, alph, bet);

if length(data) < n+1

  data = fidistinv_jacobi_setup(n, alph, bet, data);
  save_fjacobi(data, alph, bet);

end

x = zeros(size(u));
% Use data to generate an initial guess
[~,j] = histc(u, data{n+1}.u);
if any(j==0)
  error('Input values u must be between 0 and 1');
end

intervals = zeros([numel(u) 2]);
intervals(:,1) = data{n+1}.breaks(j).';

flags = (j==numel(data{n+1}.u));
intervals(~flags,2) = data{n+1}.breaks(j+1).';
intervals(flags,2) = data{n+1}.breaks(end).';

for q = 1:numel(u)
  fun = @(xx) ppval(data{n+1}, xx) - u(q);
  x(q) = fzero(fun, intervals(q,:), options);
end
