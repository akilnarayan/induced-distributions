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
N = 31;

%%%% Generate Chebyshev grid + transform
xgrid = flipud(cos(pi*(2*(0:(N-1)).'+1)/(2*N))); % equals gauss qaudrature nodes
[aa,bb] = jacobi_recurrence(N+1, -1/2, -1/2);
chebxform = poly_eval(aa, bb, xgrid, N-1);
chebxform = (diag(1./sum(chebxform.^2, 2))*chebxform).';
%%%%

fprintf('One-time setup computations: Computing induced distribution data for...\n');

% Construct piecewise polynomial data
for q = 1:length(ns)

  nn = ns(q);

  fprintf('n = %d...\n', nn);

  [a,b] = hfreud_recurrence(nn+1, alph, rho);
  xg = gauss_quadrature(a,b,nn);
  ug = idist_hfreud(xg, nn, alph, rho, 100); % Make it very accurate

  us = [0;  ug;  1 - 5*eps];
  midpts = 1/2*( us(1:end-1) + us(2:end) );
  temp = idistinv_hfreud(midpts, n, alph, rho);

  [us, inds] = sort([0; ug; 1; midpts]);
  xs = [0; xg; idistinv_hfreud(1-5*eps, n, alph, rho); temp];
  xs = xs(inds);

  dat = zeros([4+N 2*(nn+1)]);
  % Data in each col: 
  % ucenter, xcenter, exponent, scale (1-4)
  % coeffs (5-end)

  for j=1:(nn+1)
    i1 = 2*j-1;
    i2 = 2*j;

    exponents = [2/3 2/3];

    uleft = us(i1);
    ucenter = us(i2);
    uright = us(i2+1);

    xleft = xs(i1);
    xcenter = xs(i2);
    xright = xs(i2+1);

    if j == 1
      exponents(1) = rho/(rho+1);
    end
    if j==(nn+1)
      exponents(2) = 1;
    end

    ucenters(1) = uleft;
    ucenters(2) = uright;

    xcenters(1) = xleft;
    xcenters(2) = xright;

    scales(1) = (xcenter - xcenters(1)) .* abs( ucenter - ucenters(1)).^exponents(1);

    if j < nn+1
      scales(2) = (xcenters(2) - xcenter ) .* abs( ucenter - ucenters(2)).^exponents(2);
    else
      scales(2) = abs(xcenter./log(1 - ucenter).^(1/alph).*abs(1-ucenter).^exponents(2));
    end


    % i1 coeffs
    ugrid = (xgrid+1)/2*(ucenter-uleft) + uleft;
    tgrid = idistinv_hfreud(ugrid, nn, alph, rho);
    tgrid = (tgrid - xcenters(1)) .* abs(ugrid - ucenters(1)).^exponents(1);
    tgrid = tgrid/scales(1)*2 - 1;
    coeffs = chebxform*tgrid;

    dat(:,i1) = [ucenters(1); xcenters(1); exponents(1); scales(1); coeffs];

    % i2 coeffs
    ugrid = (xgrid+1)/2*(uright-ucenter) + ucenter;
    tgrid = idistinv_hfreud(ugrid, nn, alph, rho);
    if j < nn+1
      tgrid = (tgrid - xcenters(2)) .* abs(ugrid - ucenters(2)).^exponents(2);
      tgrid = tgrid/scales(2)*2 +1;
    else
      tgrid = tgrid./log(1-ugrid).^(1/alph).* abs(1 - ugrid).^exponents(2);;
      tgrid = tgrid/scales(2)*2 +1;
    end
    coeffs = chebxform*tgrid;

    dat(:,i2) = [ucenters(2); xcenters(2); exponents(2); scales(2); coeffs];

  end

  data{nn+1} = dat;

end
fprintf('Done\n');
