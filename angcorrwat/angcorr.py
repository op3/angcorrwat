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
# along with pyangdist.  If not, see <http://www.gnu.org/licenses/>.

"""Calculate angular correlations (NOT WORKING?)"""

from functools import lru_cache

from sympy import sqrt, symbols, summation, Piecewise, pi, Ynm, Eq, Mod, Abs
from sympy.physics.quantum.cg import Wigner3j

from angdist import a, f, safe_divide, flatten

@lru_cache(maxsize=None)
def isst(nu, q, l, i_n, i, sigma):
    return Piecewise(
            (0, Eq(Mod(nu, 2), 1)),
            (f(l, l, i_n, i, nu), Eq(q, 0)),
            ((-1/2 * (-1)**sigma * f(l, l, i_n, i, nu) * 
                safe_divide(Wigner3j(l, 1, l, 1, nu, -2),
                    Wigner3j(l, 1, l, -1, nu, 0))), Eq(Abs(q), 2)),
            (0, True))


def W2phot(theta1, phi1, theta2, phi2, i_i, i, i_f, e0m1, l1, e0m1p, lp1, delta1,
           l2, lp2, delta2, i_intn, lint, lintp, deltainter):
    q0, q1, q2 = symbols("q₀ q₁ q₂", integer=True)
    s = [
        [
            [
                summation(
                    ((-1)**(lam1 + lam2)/(4 * pi * sqrt(2 * lam2 + 1)) *
                    isst(lam0, q0, l1, i_i, i, e0m1) *
                    a(lint, lintp, i_intn, i, deltainter, lam1, lam2, lam0) *
                    a(l2, lp2, i_f, i_intn, delta2, lam2) *
                    Wigner3j(lam2, q2, lam1, q1, lam0, q0)) *
                    Ynm(lam1, q1, theta1, phi1) *
                    Ynm(lam2, q2, theta2, phi2),
                    (q2, -lam2, lam2),
                    (q1, -lam1, lam1),
                    (q0, -lam0, lam0))
                for lam2 in list(range(0, 5, 2))
            ] for lam1 in list(range(0, 5, 2))
        ] for lam0 in list(range(0, 5, 2))
    ]
    return sum(flatten(s))
