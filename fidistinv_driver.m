function[x] = fidistinv_driver(u, n, data)
% x = fidistinv_driver(u, n, data)
%
% Driver function for performing fast inverse distribution sampling. Uses given
% data to perform the u ---> x map associated with the order-n induced distribution. Precisely:
%
%   x(j) = F_{n(j)}^{-1} (u(j)),        1 \leq j \leq length(u)

usize = size(u);
u = reshape(u, [numel(u) 1]);
if length(u) == 0
  x = zeros(0);
  return
end

nscalar = false;
if (length(u) ~= length(n)) && (length(n) ~= 1)
  error('Inputs u and n must be the same size, or n must be a scalar')
else
  nscalar = true;
end
n = reshape(n, [numel(n) 1]);

N = max(n);
assert(length(data) >= N+1, 'Input data does not cover range of n');

x = zeros(size(u));

if length(n) > 1
  for qn = 0:N
    nmask = (n==qn);
    x(nmask) = driver_helper(u(nmask), data{qn+1});
  end
else
  x = driver_helper(u, data{n+1});
end

x = reshape(x, usize);

end

%%%%%%%%%%%%%%%% end main function

function[x] = driver_helper(u, data)
% Helper function
% u: vector
% data: cell array

tol = 1e-12; % machine eps guarding

% Construct Vandermonde-like matrix
M = size(data,1) - 6;
[aT,bT] = jacobi_recurrence(M+1, -1/2, -1/2);

edges = [-Inf data(1,:) data(2,end) Inf];
[~,j] = histc(u, edges);
B = length(edges)-1;

x = zeros(size(u));

% ``Boundary" bins: 
% 1: x ---> left edge
x(j==1) = data(3,1);
% end: x ---> right edge
x(j==B) = data(4,end);

% Interior bins:
for qb = 2:(B-1)

  umask = (j==qb);
  if not(any(umask))
    continue
  end

  q = qb-1;
  vgrid = (u(umask) - data(1,q))./(data(2,q) - data(1,q)) * 2 - 1;
  V = poly_eval(aT,bT, vgrid, M-1);

  temp = V*data(7:end,q);
  temp = temp ./ ( (1 + vgrid).^(data(5,q)) .* (1 - vgrid).^(data(6,q)) );

  if data(5,q) ~= 0
    % Set LHS to be 0
    flags = abs(u(umask) - data(1,q)) < tol;

    temp(flags) = 0;
    temp = temp * (data(4,q) - data(3,q)) + data(3,q);
  else
    % Set RHS to be 0
    flags = abs(u(umask) - data(2,q)) < tol;
    temp(flags) = 0;

    temp = temp * (data(4,q) - data(3,q)) + data(4,q);
  end

  x(umask) = temp;

end

end
