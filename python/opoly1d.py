import numpy as np
from scipy import special as sp

def jacobi_recurrence(N, alpha=0., beta=0., probability=False):
    # Returns the first N recurrence coefficient pairs for the Jacobi
    # polynomial family.

    if N < 1:
        return np.ones((0,2))

    ab = np.ones((N,2)) * np.array([beta**2.- alpha**2., 1.])

    # Special cases
    ab[0,0] = (beta - alpha) / (alpha + beta + 2.)
    ab[0,1] = np.exp( (alpha + beta + 1.) * np.log(2.) +
                      sp.gammaln(alpha + 1.) + sp.gammaln(beta + 1.) -
                      sp.gammaln(alpha + beta + 2.)
                    )

    if N > 1:
        ab[1,0] /= (2. + alpha + beta) * (4. + alpha + beta)
        ab[1,1] = 4. * (alpha + 1.) * (beta + 1.) / (
                   (alpha + beta + 2.)**2 * (alpha + beta + 3.) )

    inds = np.arange(2.,N)
    ab[2:,0] /= (2. * inds + alpha + beta) * (2 * inds + alpha + beta + 2.)
    ab[2:,1] = 4 * inds * (inds + alpha) * (inds + beta) * (inds + alpha + beta)
    ab[2:,1] /= (2. * inds + alpha + beta)**2 * (2. * inds + alpha + beta + 1.) * (2. * inds + alpha + beta - 1)

    ab[:,1] = np.sqrt(ab[:,1])

    if probability:
        ab[0,1] = 1.

    return ab

def hermite_recurrence(N, rho=0., probability=False):
    # Returns the first N recurrence coefficient pairs for the Hermite
    # polynomial family.

    if N < 1:
        return np.ones((0,2))

    ab = np.zeros((N,2))
    ab[0,1] = sp.gamma(rho+0.5)

    ab[1:,1] = 0.5*np.arange(1., N)
    ab[np.arange(N) % 2 == 1,1] += rho

    ab[:,1] = np.sqrt(ab[:,1])

    if probability:
        ab[0,1] = 1.

    return ab

def opoly1d_eval(x, n, ab):
    # Evaluates univariate orthonormal polynomials given their
    # three-term recurrence coefficients ab.

    n = np.asarray(n)
    if n.size < 1 or x.size < 1:
        return np.zeros(0)

    nmax = np.max(n)
    assert nmax < ab.shape[0]
    assert np.min(n) > -1

    p = np.zeros( x.shape + (nmax+1,) )
    xf = x.flatten()

    p[:,0] = 1/ab[0,1]

    if nmax > 0:
        p[:,1] = 1/ab[1,1] * ( (xf - ab[0,0])*p[:,0] )

    for j in range(2, nmax+1):
        p[:,j] = 1/ab[j,1] * ( (xf - ab[j,0])*p[:,j-1] - ab[j-1,1]*p[:,j-2] )

    p = p[:,n.flatten()]

    return p

if __name__ == "__main__":

    from matplotlib import pyplot as plt

    alpha = -1/2.
    beta = -1/2.
    nmax = 35
    probability_measure = True

    ab = jacobi_recurrence(nmax+1,alpha=alpha,beta=beta,probability=probability_measure)
    x = np.linspace(-1, 1, 300)

    p = opoly1d_eval(x, np.arange(nmax+1), ab)
    plt.plot(x, p[:,:10])
    plt.show()
