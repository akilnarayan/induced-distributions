function[] = save_fhfreud(data, alph, rho)
% save_fhfreud(data, alph, rho)
%
% Saves half-line Freud data to the filename given by filename_hfreud.

filename = filename_hfreud(alph, rho);
[filedir,~,~] = fileparts(filename);
if exist(filedir) ~= 7
  mkdir(filedir)
end

save(filename, 'data');
