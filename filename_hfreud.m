function[fname] = filename_hfreud(alph, rho)
% fname = filename_hfreud(alph, rho)
%
% Returns a string of value
%
%    data/half_freud_alph_bet.mat
%
% where "alph" and "rho" are replaced by their respective numerical
% values. (alph and rho are given to 4 decimal places.)

fname = fullfile('data', sprintf('half_freud_%.4f_%.4f.mat', alph, rho));
