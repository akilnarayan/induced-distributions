function[x] = fidistinv_jacobi(u, n, alph, bet);
% [x] = fidistinv_jacobi(u, n, alph, bet)
%
% A Fast Induced Distribution Inverse routine for Jacobi weights.
%
% Computes the inverse of the order-n induced primitive for the Jacobi
% distribution with parameters alph and bet. A monotone piecewise cubic
% interpolant of the distribution function is computed and saved for future
% runs.

data = load_fjacobi(n, alph, bet);

if length(data) < max(n(:))+1

  data = fidistinv_jacobi_setup(max(n(:)), alph, bet, data);
  save_fjacobi(data, alph, bet);

end

x = zeros(size(u));
%u = reshape(u, [numel(u) 1]);

if numel(n) > 1
  assert(numel(n) == numel(u));

  % Bin points with same n values
  [nsorted, inds] = sort(n(:));
  indsep = find(diff(nsorted) > 0) + 1;
else
  nsorted = n;
  inds = (1:numel(u)).';
  indsep = [];
end

% Loop over n values
for qn = 1:(length(indsep)+1)
  if qn == 1
    i1 = 1;
    if length(indsep) == 0 % This means there's only 1 n value
      i2 = numel(u);
    else
      i2 = indsep(qn)-1;
    end
  elseif qn == (length(indsep)+1)
    i1 = indsep(qn-1);
    i2 = numel(u);
  else
    i1 = indsep(qn-1);
    i2 = indsep(qn)-1;
  end
  nn = nsorted(i1);
  curr_ninds = inds(i1:i2);
  uu = u(curr_ninds);
  uu = uu(:);

  % For this n value (=nn), compute Finv(uu)
  M = size(data{nn+1},2);
  N = size(data{nn+1},1) - 4; % Number of Chebyshev coeffs

  itemp = [1, 2:2:M];
  us = data{nn+1}(1,itemp).';
  us = sort([us; 1/2*(us(1:end-1) + us(2:end))]);

  % Bin uu values into subintervals
  [~,j] = histc(uu, us);
  if any(j==0)
    error('Input values u must be between 0 and 1');
  end
  j(j==numel(us)) = numel(us) - 1;

  % Chebyshev setup + eval
  [aa,bb] = jacobi_recurrence(N+1, -1/2, -1/2);
  v = (uu - us(j))./(us(j+1) - us(j)) * 2 - 1;
  V = poly_eval(aa, bb, v, N-1);
  
  % Partition based on values of j
  [jsorted, jinds] = sort(j);
  jindsep = find(diff(jsorted) > 0) + 1;

  xx = zeros(size(uu));

  % Loop over subintervals
  for jj = 1:(length(jindsep)+1)
    if jj == 1
      i1 = 1;
      if length(jindsep) == 0
        i2 = numel(uu);
      else
        i2 = jindsep(jj)-1;
      end
    elseif jj == (length(jindsep)+1)
      i1 = jindsep(jj-1);
      i2 = numel(uu);
    else
      i1 = jindsep(jj-1);
      i2 = jindsep(jj)-1;
    end

    j = jsorted(i1);
    currinds = jinds(i1:i2);

    temp = V(currinds,:)*data{nn+1}(5:end,j);
    if mod(j,2) == 0
      flags = abs(v(currinds)-1) < 1e-8; % These are really close to uright
      sgn = -1;
    else
      flags = abs(v(currinds)+1) < 1e-8; % There are really close to uleft
      sgn = +1;
    end

    temp = (temp + sgn)/2*data{nn+1}(4,j) ./ abs( uu(currinds) - data{nn+1}(1,j)).^data{nn+1}(3,j);
    temp = temp + data{nn+1}(2,j);

    temp(flags) = data{nn+1}(2,j); % Peg to x endpoint

    xx(currinds) = temp;

  end

  x(curr_ninds) = xx;

end

% Make sure nothing is outside [-1,1]
flags = abs(x) > 1;
x(flags) = sign(x(flags));

% For Half Freud: for large x the distribution looks like x^{2n+rho} exp(-x^alph) * x^(1-alph).
