function[a,b] = quadratic_modification_C(alph, bet, x0)
% quadratic_modification_C -- Modifies recurrence coefficients
%
% [a,b] = quadratic_modification_C(alph, bet, x0)
%
%   Performs a quadratic modification of orthogonal polynomial recurrence
%   coefficients. The inputs alph and bet are three-term recurrence
%   coefficients for a d(mu)-orthonormal polynomial family p_n:
%
%     sqrt(bet_{n+1}) p_{n+1} = (x - alph_n) p_n - sqrt(bet_n) p_{n-1}
%
%   This function transforms the alph, bet into new coefficients a, b such that
%   the new coefficients determine a polynomial family that is orthonormal
%   under the weight (x-x0)^2*d(mu).
%
%   This function uses the q functions to perform the recurrence updates,
%   instead of the standard orthogonal polynomial basis.
%
%   Outputs:
%       a: column vector (size (N-2)) of recurrence coefficients
%       b: column vector (size (N-2)) of recurrence coefficients
%
%   Inputs:
%       x0   : real scalar
%       alph : column vector (length-N) of recurrence coefficients
%       bet  : column vector (length-N) of reucrrence coefficients

N = length(alph);
assert( length(bet) == N );
assert(N > 2);

alph = alph(:);
bet = bet(:);

% Output recurrence coefficients
a = zeros([N-2 1]);
b = zeros([N-2 1]);

C = reshape(C_eval(alph, bet, x0, N-1), [N 1]);

% q is length N --- new coefficients have length N-2
acorrect = zeros([N-2 1]);
bcorrect = zeros([N-2 1]);

temp1 = sqrt(bet(2:N)).*C(2:N).*C(1:(N-1))./sqrt(1 + C(1:(N-1)).^2);
temp1(1) = sqrt(bet(2))*C(2); % special case

acorrect = diff(temp1);

temp1 = 1 + C(1:(N-1)).^2;
bcorrect = temp1(2:end)./temp1(1:(end-1));
bcorrect(1) = (1 + C(2).^2)./C(1).^2;

b = bet(2:N-1) .* bcorrect;
a = alph(2:N-1) + acorrect;
