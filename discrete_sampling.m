function[x] = discrete_sampling(N, prob, states)
% discrete_sampling -- samples iid from a discrete probability measure
%
% x = discrete_sampling(N, prob, states)
%
%   Generates N iid samples from a random variable X whose probability mass
%   function is 
%
%     prob(X = states(j)) = prob(j),    1 <= j <= length(prob).
%
%   If states is not given, the states are gives by 1 <= state <= length(prob).

p = prob(:)/sum(prob(:));

[~,bins] = histc(rand(N,1), [0;cumsum(p)]);

if nargin < 3
  x = bins;
else
  assert(numel(states) == numel(prob));
  x = states(bins);
end
