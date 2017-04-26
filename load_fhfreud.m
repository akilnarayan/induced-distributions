function[data] = load_fhfreud(n, alph, rho)
% [data] = load_fhfreud(n, alph, rho)
%
% Loads available data for computing induced half-line Freud distribution inverses.
% Attempts to load this data from the "data" subdirectory, looking for a
% filename given by filename_hfreud.
%
% Returns an empty cell array if either the subdirectory "data" does not exist
% or the desired filename does not exist.

filename = filename_hfreud(alph, rho);

if exist(filename) == 2;
  data = load(filename);
  data = data.data;
else
  data = cell(0);
end
