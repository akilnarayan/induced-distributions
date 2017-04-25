function[data] = load_fjacobi(n, alph, bet)
% [data] = load_fjacobi(n, alph, bet)
%
% Loads available data for computing induced Jacobi primitive inverses.
% Attempts to load this data from the "data" subdirectory, looking for a
% filename given by jacobi_data_filename.
%
% Returns an empty cell array if either the subdirectory "data" does not exist
% or the desired filename does not exist.

filename = filename_jacobi(alph, bet);

if exist(filename) == 2;
  data = load(filename);
  data = data.data;
else
  data = cell(0);
end
