function[x] = fidistinv_hfreud(u, n, alph, rho)
% [x] = fidistinv_hfreud(u, n, alph, rho)
%
% A Fast Induced Distribution Inverse routine for Half-line Freud weights.
%
% Computes the inverse of the order-n induced primitive for the half-line Freud
% distribution with parameters alph and rho.

data = load_fhfreud(n, alph, rho);

if length(data) < max(n(:))+1

  data = fidistinv_hfreud_setup(max(n(:)), alph, rho, data);
  save_fhfreud(data, alph, rho);

end

x = zeros(size(u));
if numel(u) == 0
  return
end
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

    if j < 2*(nn+1)
      temp = (temp + sgn)/2*data{nn+1}(4,j) ./ abs( uu(currinds) - data{nn+1}(1,j)).^data{nn+1}(3,j);
      temp = temp + data{nn+1}(2,j);

    else
      temp(temp > 1) = 1-1e-5; % TODO: This is a bad fix based on a poor approximation at RHS
      temp = (temp + sgn)/2*data{nn+1}(4,j) .* log( 1 - uu(currinds) ).^(1/alph) ./ abs(1 - uu(currinds)).^data{nn+1}(3,j);
    end

    temp(flags) = data{nn+1}(2,j); % Peg to x endpoint
    xx(currinds) = temp;

  end

  x(curr_ninds) = xx;

end

% Make sure nothing is negative
flags = x < 0;
x(flags) = 0;

% For Half Freud: for large x the distribution looks like x^{2n+rho} exp(-x^alph) * x^(1-alph).
