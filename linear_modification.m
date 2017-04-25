function[a,b] = linear_modification(alph, bet, x0)
% linear_modification -- Modifies recurrence coefficients
%
% [a,b] = linear_modification(alph, bet, x0)
%
%   Performs a linear modification of orthogonal polynomial recurrence
%   coefficients. The inputs alph and bet are three-term recurrence
%   coefficients for a d(mu)-orthonormal polynomial family p_n:
%
%     sqrt(bet_{n+1}) p_{n+1} = (x - alph_n) p_n - sqrt(bet_n) p_{n-1}
%
%   This function transforms the alph, bet into new coefficients a, b such that
%   the new coefficients determine a polynomial family that is orthonormal
%   under the weight |x-x0|*d(mu), when x0 \not\in \supp \mu.
%
%   The appropriate sign of the modification (+/- (x-x0)) is inferred from the
%   sign of (alph(1) - x0). Since alph(1) is the zero of p_1, then it is in
%   \supp \mu.
%
%   Outputs:
%       a: column vector (size (N-1)) of recurrence coefficients
%       b: column vector (size (N-1)) of recurrence coefficients
%
%   Inputs:
%       x0   : real scalar, assumed \not \in \supp \mu
%       alph : column vector (length-N) of recurrence coefficients
%       bet  : column vector (length-N) of reucrrence coefficients
%       x0   : real scalar, +1 or -1, so that sgn*(x-x0) > 0 on \supp \mu

N = length(alph);
assert( length(bet) == N );
assert(N > 1);

sgn = sign(alph(1) - x0);

r = reshape(abs(ratio_eval(alph, bet, x0, N-1)), size(alph(1:end-1)));
% r is length N-1

acorrect = zeros([N-1 1]);
bcorrect = zeros([N-1 1]);

acorrect = sqrt(bet(2:N))./r;
acorrect(2:end) = diff(acorrect);

bcorrect = sqrt(bet(2:N)).*r;
bcorrect(2:end) = bcorrect(2:end)./bcorrect(1:end-1);

b = bet(1:N-1) .* bcorrect(1:N-1);
a = alph(1:N-1) + sgn*acorrect(1:N-1);
