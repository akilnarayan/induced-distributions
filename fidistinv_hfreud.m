function[x] = fidistinv_hfreud(u, n, alph, rho)
% [x] = fidistinv_hfreud(u, n, alph, rho)
%
% A Fast Induced Distribution Inverse routine for half-Freud weights.
%
% Computes the inverse of the order-n induced primitive for the half-Freud
% distribution with parameters alph and rho. A monotone piecewise cubic
% interpolant of the distribution function is computed and saved for future
% runs.

persistent options
if isempty(options)
  options = optimset(@fzero);
  options.Display = 'off';
  options.TolFun = 10*eps;
  options.TolX = 10*eps;
end

data = load_fhfreud(n, alph, rho);

if length(data) < n+1

  data = fidistinv_hfreud_setup(n, alph, rho, data);
  save_fhfreud(data, alph, rho);

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
