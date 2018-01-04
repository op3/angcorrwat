"""Calculate angular distributions"""

from functools import lru_cache

from sympy import sqrt, factorial, cos, Symbol, summation, Piecewise, pi, prod
from sympy.physics.quantum.cg import CG, Wigner3j, Wigner6j, Wigner9j
from sympy.functions import legendre, assoc_legendre

deg = pi/180.


def legenp(n, m, x):
    return Piecewise((0, m > n), (assoc_legendre(n, m, x), True))


def safe_divide(numerator, denominator):
    return denominator and numerator/denominator


def flatten(items, seqtypes=(list, tuple)):
    for i, _ in enumerate(items):
        while i < len(items) and isinstance(items[i], seqtypes):
            items[i:i+1] = items[i]
    return items


def kappa(nu, l, lp):
    return Piecewise(((-sqrt(factorial(nu - 2)/factorial(nu + 2)) *
                       CG(l, 1, lp, 1, nu, 2) /
                       CG(l, 1, lp, -1, nu, 0)),
                      ((nu > 1) & (nu < l + lp + 1))), (0, True))


@lru_cache(maxsize=None)
def f(l, lp, i2, i1, e, f=None, g=None):
    if f is None:
        return _f_5(l, lp, i2, i1, e)
    elif g is None:
        return _f_6(l, lp, i2, i1, e, f)
    else:
        return _f_gen(l, lp, i2, i1, e, f, g)


def _f_5(l, lp, i2, i1, nu):
    return Piecewise((
        ((-1)**(i2 + i1 - 1) *
         sqrt((2 * l + 1) * (2 * lp + 1) * (2 * i1 + 1) * (2 * nu + 1)) *
         Wigner3j(l, 1, lp, -1, nu, 0) * Wigner6j(l, lp, nu, i1, i1, i2)),
        ((l + lp >= nu) & (abs(l-lp) <= nu) &
         (lp + i1 >= i2) & (abs(lp - i1) <= i2) &
         (l + i2 >= i1) & (abs(l - i2) <= i1) &
         (abs(nu - i1) <= i1))), (0, True))


def _f_6(l, lp, i2, i1, lam2, lam1):
    return Piecewise((
        ((-1)**(i1 + i2 + l) *
         sqrt((2 * i1 + 1) * (2 * i2 + 1) * (2 * lam1 + 1)) *
         Wigner6j(i1, i1, lam1, i2, i2, l)),
        ((l == lp) & (lam1 == lam2) &
         (i1 + i1 >= lam1) & (abs(i1 - i1) <= lam1) &
         (i1 + i2 >= l) & (abs(i1 - i2) <= l) &
         (i1 + l >= i2) & (abs(i1 - l) <= i2) &
         (abs(lam1 - i2) <= i2))), (0, True))


def _f_gen(l, lp, i2, i1, nu, lam2, lam1):
    if(nu > l + lp or nu > lam1 + lam2 or (lam1 > 2*i1 and lam2 > 2*i2)):
        return 0
    return (sqrt(prod([2*i + 1 for i in locals().values()])) *
            (-1)**(lp + lam1 + lam2 + 1) *
            Wigner3j(l, 1, lp, -1, nu, 0) * 
            Wigner9j(i2, l, i1, i2, lp, i1, lam2, nu, lam1))


def u(lam1, l, lp, i2, i1, deltainter):
    return ((f(l, l, i2, i1, lam1, lam1) +
             2 * deltainter * f(l, lp, i2, i1, lam1, lam1) +
             deltainter**2 * f(lp, lp, i2, i1, lam1, lam1)) /
            ((1 + deltainter**2) * sqrt(2 * lam1 + 1)))


def a(l, lp, i_n, i, delta, nu, lam2=None, lam1=None):
    if lam2 is None and lam1 is None:
        return ((f(l, l, i_n, i, nu) + 2 * delta * f(l, lp, i_n, i, nu) +
                 delta**2 * f(lp, lp, i_n, i, nu))/(1 + delta**2))
    return ((f(l, l, i_n, i, nu, lam2, lam1) +
             2 * delta * f(l, lp, i_n, i, nu, lam2, lam1) +
             delta**2 * f(lp, lp, i_n, i, nu, lam2, lam1))/(1 + delta**2))


def b(l, lp, i_n, i, delta, nu):
    return ((f(l, l, i_n, i, nu) + (-1)**(l + lp) * 2 *
             delta * f(l, lp, i_n, i, nu) + delta**2 *
             f(lp, lp, i_n, i, nu)) /
            (1 + delta**2))


def bp(nu, theta, phi, sigma, l, sigmap, lp, i_n, i, delta):
    return (b(l, lp, i_n, i, delta, nu) * legendre(nu, cos(theta)) +
            1/(1 + delta**2) * cos(2 * phi) * legenp(nu, 2, cos(theta)) *
            (f(l, l, i_n, i, nu) * (-1)**sigma * kappa(nu, l, l) +
             (-1)**(l + lp) * 2 * delta * f(l, lp, i_n, i, nu) *
             (-1)**sigmap * kappa(nu, l, lp) +
             delta**2 * f(lp, lp, i_n, i, nu) *
             (-1)**sigmap * kappa(nu, lp, lp)))


@lru_cache(maxsize=None)
def isst(nu, q, l, i_n, i, sigma):
    if nu % 2:
        return 0
    elif q == 0:
        return f(l, l, i_n, i, nu)
    elif abs(q) == 2:
        return (-1/2 * (-1)**sigma * f(l, l, i_n, i, nu) *
                safe_divide(Wigner3j(l, 1, l, 1, nu, -2),
                            Wigner3j(l, 1, l, -1, nu, 0)))
    else:
        return 0


"""
def _W(theta, phi,
       i_i, sigma, l1, sigmap, lp1, delta1,
       [i2, l2, lp2, delta2],
       i_f):
"""
def _W(theta, phi,
       i_i, i, i_f, sigma, l1, sigmap, lp1, delta1, l2, lp2, delta2,
       iintn=None, lint=None, lintp=None, deltaint=None):
    """
    Spin sequence is: i_i (sigma l1, sigmap lp1) i (l2, lp2) iintn (lint, lintp) i_f
    """
    if iintn is not None:
        raise NotImplementedError()
    nu = Symbol('nu', integer=True)
    return summation(
        bp(2 * nu, theta, phi, sigma, l1, sigmap, lp1, i_i, i, delta1) *
        a(l2, lp2, i_f, i, delta2, 2 * nu), (nu, 0, 2))


def W(*args, **kwargs):
    """
    Angular distribution for two successive γ-rays.
    The incident γ-ray is assumed to be polarized. For the calculation of
    unpolarized incident γ-rays, horizontal (φ = 0) and vertical (φ = π/2)
    polarizations must be added up.

    Spin sequence is: i_i (sigma l1, sigmap lp1) i (l2, lp2) i_f
    For the first gamma ray due to the polarization assumption the electric or
    magnetic character of the radiation must be known.

    sigma_i l_i is the multipole radiation of the transition mixing  with
    sigma_i l_i' (e.g., M1/E2 mixing) with the multipole mixing ratio delta_i.

    Args:
        theta: angle to z-axis defined by first γ
        phi: angle rotating around z-axis
        i_i: spin of initial state
        i : spin of intermediate state
        i_f: spin of final state.
        sigma_i: electric (0) or magnetic (1) character of incident γ.
        l_i: multipole order of the radiation
        delta_i: multipole mixing ratio (convention of Krane/Steffen/Wheeler)
    """
    return _W(*args, **kwargs).doit().simplify()


def analyzing_power(*args, **kwargs):
    """
    The analyzing power defined as (W_hor - W_ver)/(W_hor + W_ver)
    """
    return ((W(pi/2, 0, *args, **kwargs) - W(pi/2, pi/2, *args, **kwargs)) /
            (W(pi/2, 0, *args, **kwargs) + W(pi/2, pi/2, *args, **kwargs))
            ).doit().simplify()
