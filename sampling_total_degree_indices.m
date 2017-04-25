function[lambdas] = sampling_total_degree_indices(N, d, k)
% sampling_total_degree_indices -- generates random multi-indices
%
% lambdas = sampling_total_degree_indices(N, d, k)
%
%   Chooses N random multi-indices (with the uniform probability law) from the
%   set of d-variate multi-indices whose total degree is k and less.
%
%   The output lambdas is an N x d matrix, with each row containing one of
%   these multi-indices.
%
%   This function draws N samples without actually generating the full total
%   degree set.

lambdas = zeros([N d]);

degrees = discrete_sampling(N, pdjk(d, k), 0:k).';

for i = 1:(d-1)
  for n = 1:N
    lambdas(n,i) = discrete_sampling(1, pdjk(d-i, degrees(n)), degrees(n):-1:0);
  end

  degrees = degrees - lambdas(:,i);
end

lambdas(:,d) = degrees;

end

function[p] = pdjk(d,k)

  j = 0:k;
  p = exp( log(d) + gammaln(k+1) - gammaln(j+1) + gammaln(j+d) - gammaln(k + d + 1) );

  assert(abs(sum(p) - 1) < 1e-8);

end
