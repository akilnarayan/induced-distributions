function[data] = fidistinv_hfreud_setup(n, alph, rho, data)
% data = fidistinv_hfreud_setup(n, alph, rho, data)
%
% Computes coefficients for a piecewise cubic interpolant of the degree-n
% induced distribution for half-Freud weights. These coefficients are appended
% to the data cell array, inside data{n+1}.
%
% This function also computes the "downward closed" data, that is, if
% length(data)==k, meaning data up to degree k-1 is present, then this function
% computes and stores data for degrees k, k+1, ..., n.

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

fprintf('One-time setup computations: Computing primitive data for...\n');

% Construct piecewise polynomial data
for q = 1:length(ns)

  nn = ns(q);

  fprintf('n = %d...\n', nn);
  
  % The data points we'll choose are as follows:
  % 2*nn+1 equally spaced data on u
  % Q data points on each interval bounded by zeros of p_nn (Q*(nn+1) + 1 points)

  [a,b] = hfreud_recurrence(nn+1, alph, rho);
  x = gauss_quadrature(a, b, nn);


  x = [0; x; (maxapprox_hfreud(alph, rho, nn):hfreud_tolerance(nn, alph, rho, eps)).'];
  %x = [0; x; (half_freud_mrs_max_guess(alph, rho, n):half_freud_effective_limit(n, alph, rho, 1e-16)).'];

  % Create Q-1 points per interval
  xn = diff(x) * xtemplate + repmat(x(1:end-1), [1 Q-1]);
  xn = [reshape(xn.', [numel(xn) 1]); x(end)];
  X = numel(xn);

  un = zeros(size(xn));
  un = idist_hfreud(xn, nn, alph, rho);
  %for qq = 1:X
  %  un(qq) = gfreud_induced_primitive(xn(qq), nn, alph, rho);
  %end

  % Now nn*R data points equally-spaced in u space
  temp = linspace(0, 1, R + 2).';
  temp([1 end]) = [];

  un = [un; temp];
  %xn = [xn; gfreud_primitive_inverse(temp, nn, alph, rho)];
  xn = [xn; idistinv_hfreud(temp, nn, alph, rho)];

  % Now nn*R data points equally-spaced in x space
  temp = linspace(x(1), x(end), R+2).';
  temp([1 end]) = [];

  X = numel(xn);
  xn = [xn; temp];
  un = [un; zeros(size(temp))];
  un((X+1):numel(xn)) = idist_hfreud(xn((X+1):numel(xn)), nn, alph, rho);
  %for qq = (X+1):numel(xn);
  %  un(qq) = gfreud_induced_primitive(xn(qq), nn, alph, rho);
  %end

  % Sort
  [xn, inds] = sort(xn);
  un = un(inds);

  % Some non-monotonic behavior happens very close to u=1. Remove these.
  while true
    inds = [false; diff(un) < 0];
    if not(any(inds))
      break
    end
    un(inds) = [];
    xn(inds) = [];
  end

  % Remove identical x values
  while true
    inds = [false; diff(xn) <= 0];
    if not(any(inds))
      break
    end
    un(inds) = [];
    xn(inds) = [];
  end

  data{nn+1} = pchip(xn, un);
  data{nn+1}.u = un;

end
fprintf('Done\n');
