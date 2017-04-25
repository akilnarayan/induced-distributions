function[fname] = filename_jacobi(alph, bet)
% fname = filename_jacobi(alph, bet)
%
% Returns a string of value
%
%    data/jacobi_alph_bet.mat
%
% where "alph" and "bet" are replaced by their respective numerical
% values. (alph and bet are given to 4 decimal places.)

fname = fullfile('data', sprintf('jacobi_%.4f_%.4f.mat', alph, bet));
