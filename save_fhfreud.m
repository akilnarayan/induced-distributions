function[] = save_fhfreud(data, alph, rho)
% save_fhfreud(data, alph, rho)
%
% Saves half-line Freud data to the filename given by filename_hfreud.

save(filename_hfreud(alph, rho), 'data');
