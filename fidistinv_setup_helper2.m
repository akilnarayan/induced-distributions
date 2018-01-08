function[data] = fidistinv_setup_helper2(ug, idistinv, exponents, M);
% [data] = fidistinv_setup_helper2(ug, idistinv, exponents, iV)
%
% Helper function for setup computations for fast induced distribution
% inversion.
%
% Inputs:
%  ug: vector, values of the induced distribution
%  idistinv: function handle, with calling syntax idistinv(ug)
%  exponents: array of exponents
%  iV: size of Chebyshev transform to take

% Chebyshev transform setup
vgrid = cos(linspace(pi, 0, M)).';
[aT,bT] = jacobi_recurrence(M+1, -1/2, -1/2);
V = poly_eval(aT, bT, vgrid, M-1);
iV = inv(V);

ugrid = zeros([M length(ug)-1]);
xgrid = zeros([M length(ug)-1]);
xcoeffs = zeros([M length(ug)-1]);
for q = 1:(length(ug)-1)

  ugrid(:,q) = (vgrid + 1)/2 * (ug(q+1) - ug(q)) + ug(q);
  if q == length(ug)-1
    xgrid(:,q) = idistinv(ugrid(:,q));
  else
    xgrid(:,q) = idistinv(ugrid(:,q));
  end

  temp = xgrid(:,q);
  if exponents(1,q) ~= 0
    temp = (temp - xgrid(1,q))/(xgrid(end,q) - xgrid(1,q));
  else
    temp = (temp - xgrid(end,q))/(xgrid(end,q) - xgrid(1,q));
  end

  temp = temp .* (1 + vgrid).^(exponents(1,q)) .* (1 - vgrid).^(exponents(2,q));

  % The assumption below is that we choose exponents so that the function is
  % linear and has a zero at endpoints.
  temp(~isfinite(temp)) = 0;

  xcoeffs(:,q) = iV*temp;

end

data = zeros([M + 6, length(ug)-1]);
for q = 1:(length(ug)-1)
  data(:,q) = [ug(q); ug(q+1); xgrid(1,q); xgrid(end,q); exponents(:,q); xcoeffs(:,q)];
end

