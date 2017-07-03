function[x] = idist_mixture_sampling(varargin)
% x = idist_mixture_sampling(M, Lambdas, univ_inv)
% x = idist_mixture_sampling(Lambdas, univ_inv)
%
% Performs tensorial inverse transform sampling from an additive mixture of
% tensorial induced distributions, generating M samples. The measure this
% samples from is the order-Lambdas induced measure, which is an additive
% mixture of tensorial measures. Each tensorial measure is defined a row of
% Lambdas. The dimension d of the problem is inferred from size(Lambdas, 2).
%
% The function univ_inv is a function handle that inverts a univariate order-n
% induced distribution. It must support the syntax
%
%    univ_inv( u, n )
%
% where u can be any size array of elements between 0 and 1, and n is a
% non-negative integer.
%
% The second calling syntax sets M = size(Lambdas, 1), so that Lambdas already
% repesents randomly-generated multi-indices.

if nargin < 3
  Lambdas = varargin{1};
  univ_inv = varargin{2};

  [M, d] = size(Lambdas);
else
  M = varargin{1};
  Lambdas = varargin{2};
  univ_inv = varargin{3};

  [K, d] = size(Lambdas);
end

assert(M > 0);

x = zeros([M d]);

if nargin < 3
  %degrees = Lambdas(:);
else
  % This selects M random samples on 1, 2, ..., K
  ks = ceil( K*rand([M 1]) );
  ks(ks > K) = K;
  Lambdas = Lambdas(ks,:);

  %degrees = reshape(Lambdas(ks, :), [M*d 1]);
end

x = univ_inv(rand([M d]), Lambdas);

%for j = 1:d
%
%  x(:,j) = univ_inv(rand([M 1]), Lambdas(:,j));
%
%end

% Basic strategy:
%   - tabulate all univariate degrees
%   - bin degrees
%   - for each degree n, generate order-n idist samples
%   - place those samples in appropriate locations

% We need a max degree count so that we can do univariate sampling
%max_degree = max(Lambdas(:));
%
%[degree_sample_count, bin_indices] = histc( degrees , (0:(max_degree+1)) - 0.5 );
%assert(degree_sample_count(end) == 0);
%assert(sum(degree_sample_count) == M*d); % Number of samples necessary
%degree_sample_count(end) = [];
%
%for n = 0:max_degree
%  u = rand([degree_sample_count(n+1) 1]);
%  if numel(u) > 0
%    inds = find(bin_indices==n+1);
%
%    % Sample and place in appropriate locations in x
%    x(inds) = univ_inv(u, n);
%  end
%end
