function[] = save_fjacobi(data, alph, bet)
% save_fjacobi(data, alph, bet)
%
% Saves Jacobi data to the filename given by jacobi_data_filename.

filename = filename_jacobi(alph, bet);
[filedir,~,~] = fileparts(filename);
if exist(filedir) ~= 7
  mkdir(filedir)
end

save(filename, 'data');
