function[lambdas] = total_degree_indices(d, k)
% lambdas = total_degree_indices(d, k)
%
% Generates all d-variate multi-indices whose entries sum to k or less. The
% total number of these indices is N = nchoosek(d+k, k).
%
% The output lambdas is an N x d array, where each row contains a multi-index.
% The indices are ordered via graded reverse lexicographic ordering, where the
% grade is total degree.

lambdas = grlex(1:nchoosek(d+k, d), d);
