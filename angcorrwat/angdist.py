#  SPDX-License-Identifier: GPL-3.0+
#
# Copyright © 2017-2018 O. Papst.
#
# This file is part of angcorrwat.
#
# angcorrwat is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# angcorrwat is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with angcorrwat.  If not, see <http://www.gnu.org/licenses/>.

"""Calculate angular distributions"""

from functools import lru_cache
import inspect

from sympy import sqrt, factorial, cos, Symbol, summation, Piecewise, pi, prod
from sympy.physics.quantum.cg import CG, Wigner3j, Wigner6j, Wigner9j
from sympy.functions import legendre, assoc_legendre


def wrap_function(function):
    """Create a function that copies another functions docs and signature"""
    def _(func):
        func.__doc__ = function.__doc__
        func.__signature__ = inspect.signature(function)
        return func
    return _


def safe_divide(numerator, denominator):
    """Division by 0 returns 0"""
    return denominator and numerator/denominator


def flatten(items, seqtypes=(list, tuple)):
    """Turn nested lists into a single flat list"""
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
    if g is None:
        return _f_6(l, lp, i2, i1, e, f)
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
    return Piecewise(
        (0, (nu > l + lp) | (nu > lam1 + lam2) | ((lam1 > 2 * i2) & (lam2 > 2*i2))),
        ((
            sqrt(prod([2*i + 1 for i in locals().values()])) *
            (-1)**(lp + lam1 + lam2 + 1) *
            Wigner3j(l, 1, lp, -1, nu, 0) *
            Wigner9j(i2, l, i1, i2, lp, i1, lam2, nu, lam1)), True))


def u(lam1, l, lp, i2, i1, deltainter):
    return ((f(l, l, i2, i1, lam1, lam1) +
             2 * deltainter * f(l, lp, i2, i1, lam1, lam1) +
             deltainter**2 * f(lp, lp, i2, i1, lam1, lam1)) /
            ((1 + deltainter**2) * sqrt(2 * lam1 + 1)))


@lru_cache(maxsize=None)
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
            1/(1 + delta**2) * cos(2 * phi) * assoc_legendre(nu, 2, cos(theta)) *
            (f(l, l, i_n, i, nu) * (-1)**sigma * kappa(nu, l, l) +
             (-1)**(l + lp) * 2 * delta * f(l, lp, i_n, i, nu) *
             (-1)**sigmap * kappa(nu, l, lp) +
             delta**2 * f(lp, lp, i_n, i, nu) *
             (-1)**sigmap * kappa(nu, lp, lp)))


def _W(theta, phi, initial_state, excited_state, cascade):
    """
    Angular distribution of γ-radiation emitted by an excited nucleus.
    The incident γ-ray is assumed to be polarized. For the calculation of
    unpolarized incident γ-rays, horizontal (φ = 0) and vertical (φ = π/2)
    polarizations must be added up.

    Args:
        theta: angle to z-axis defined by first γ
        phi: angle rotating around z-axis
        initial_state: state that is initially populated before excitation
            (normally the ground state). Format is
            [total angular momentum (J), parity (π)].
        excited_state: state that is excited by the γ-ray. Format is
            [total angular momentum (J), parity (π),
                multipole mixing ratio (δ)].
        cascade: List of states that occur in the subsequent cascade. Format is
            [[total angular momentum (J), multipole mixing ratio (δ)], …].
            For the multipole mixing ratio, the convention by
            Krane/Steffen/Wheeler is used. Only the last state is observed.

    Returns:
        The angular distribution of the given cascade. Might still
        depend on input variables.
    """
    # 0: initial state
    # ex: excited state
    # prev: previous state
    # int: intermediate states
    # final: final state
    nu = Symbol('nu', integer=True)

    tot_ang_mom_0, parity_0 = initial_state
    tot_ang_mom_ex, parity_ex, delta_ex = excited_state
    tot_ang_mom_final, delta_final = cascade[-1]

    orbit_ang_mom_ex = max(abs(tot_ang_mom_ex - tot_ang_mom_0), 1)
    orbit_ang_mom_exp = orbit_ang_mom_ex + 1
    orbit_ang_mom_ex_sigma = (orbit_ang_mom_ex + parity_0 + parity_ex) % 2
    orbit_ang_mom_exp_sigma = (orbit_ang_mom_exp + parity_0 + parity_ex) % 2

    middle = 1
    tot_ang_mom_prev = tot_ang_mom_ex
    for state in cascade[:-1]:
        tot_ang_mom_int, delta_int = state
        orbit_ang_mom_int = max(abs(tot_ang_mom_prev - tot_ang_mom_int), 1)
        orbit_ang_mom_intp = orbit_ang_mom_int + 1

        middle *= u(2 * nu, orbit_ang_mom_int, orbit_ang_mom_intp,
                    tot_ang_mom_int, tot_ang_mom_prev, delta_int)
        tot_ang_mom_prev = tot_ang_mom_int

    orbit_ang_mom_final = max(abs(tot_ang_mom_prev - tot_ang_mom_final), 1)
    orbit_ang_mom_finalp = orbit_ang_mom_final + 1

    return summation(
        bp(2 * nu, theta, phi, orbit_ang_mom_ex_sigma, orbit_ang_mom_ex,
           orbit_ang_mom_exp_sigma, orbit_ang_mom_exp, tot_ang_mom_0,
           tot_ang_mom_ex, delta_ex) *
        middle *
        a(orbit_ang_mom_final, orbit_ang_mom_finalp, tot_ang_mom_final,
          tot_ang_mom_prev, delta_final, 2 * nu),
        (nu, 0, 2))


@wrap_function(_W)
def W(*args, **kwargs):
    return _W(*args, **kwargs).doit().simplify()


def analyzing_power(*args, **kwargs):
    """
    The analyzing power defined as (W_hor - W_ver)/(W_hor + W_ver)
    """
    return ((W(pi/2, 0, *args, **kwargs) - W(pi/2, pi/2, *args, **kwargs)) /
            (W(pi/2, 0, *args, **kwargs) + W(pi/2, pi/2, *args, **kwargs))
            ).doit().simplify()
