function[] = save_fjacobi(data, alph, bet)
% save_fjacobi(data, alph, bet)
%
% Saves Jacobi data to the filename given by jacobi_data_filename.

save(filename_jacobi(alph, bet), 'data');
