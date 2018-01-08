function[ug, exponents] = fidistinv_setup_helper1(ug, exps)
% [ug, exponents] = fidistinv_setup_helper1(ug, exps)
%
% Helper function for setup computations for fast induced distribution
% inversion.
%
% Inputs:
%  ug: vector, values of the induced distribution at the Gauss quadrature
%      nodes, plus endpoint values (usually 0 and 1).
%  exps: a size-2 vector containing custom exponent values at the endpoints 
%        I.e., at (u=0 and 1, respectively.)

% Sample midpoints
ug_mid = 1/2 * ( ug(1:(end-1)) + ug(2:end) );
ug = sort([ug; ug_mid]);

exponents = zeros([2, size(ug,1)-1]);
for q = 1:(size(ug,1)-1)

  if mod(q,2) == 1
    exponents(1,q) = 2/3;
  else
    exponents(2,q) = 2/3;
  end

end

% Exponents for endpoints
exponents(1,1) = exps(1);
exponents(2,end) = exps(2);
