function[F] = hermite_equilibrium_distribution_conjecture(u, d)
% hermite_equilibrium_distribution -- Radial distribution function for Hermite polynomials
%
% F = hermite_equilibrium_distribution_conjecture(u, d)
%
%   Computes the distribution function F(u) for u \in [0,1] for the
%   (sqrt(2*degree)-scaled) radius coordinate of the d-dimension equilibrium
%   distrubtion.
%
%   This computation is based on a conjecture for the density of the radius,
%   which is
%
%      f_r(u) = sqrt(1 - u^2)^(d/2)

F = betainc(u.^2, d/2, 1 + d/2);
