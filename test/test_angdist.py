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

import pytest

from sympy import Symbol, cos, sin
from angcorrwat import W

theta = Symbol('theta')
phi = Symbol('phi')
delta = Symbol('delta')


def is_equal(expression1, expression2):
    return (expression1 - expression2).simplify() == 0


@pytest.mark.parametrize("description, W_args, result", [
    ["0+ → 1- → 0+", [theta, phi, [0, 1], [1, 0, 0], [[0, 0]]],
     ((1 + cos(theta)**2 - cos(2 * phi) * sin(theta)**2) * 3 / 4)],
    ["0+ → 1+ → 0+", [theta, phi, [0, 1], [1, 1, 0], [[0, 0]]],
     ((1 + cos(theta)**2 + cos(2 * phi) * sin(theta)**2) * 3 / 4)],
    ["0+ → 2+ → 0+", [theta, phi, [0, 1], [2, 1, 0], [[0, 0]]],
     ((2 + cos(2 * theta) + cos(4 * theta) - 2 * cos(2 * phi) *
       (1 + 2 * cos(2 * theta)) * sin(theta)**2) * 5 / 8)],
    ["0+ → 1- → 2+", [theta, phi, [0, 1], [1, 0, 0], [[2, 0]]],
     ((13 + cos(theta)**2 - cos(2 * phi) * sin(theta)**2) * 3 / 40)],
    ["0+ → 1+ → 2+", [theta, phi, [0, 1], [1, 1, 0], [[2, 0]]],
     ((13 + cos(theta)**2 + cos(2 * phi) * sin(theta)**2) * 3 / 40)],
    ["0+ → 2+ → 2+", [theta, phi, [0, 1], [2, 1, 0], [[2, 0]]],
     ((7 + 3 * cos(theta)**2 + 3 * cos(2 * phi) * sin(theta)**2) * 1 / 8)],
    ])
def test_single(description, W_args, result):
    print(f"Calculating angular distribution {description}.")
    assert is_equal(W(*W_args), result)


@pytest.mark.parametrize("description, W_args, result", [
    ["0+ → 1- → 2+ → 0+", [theta, phi, [0, 1], [1, 0, 0], [[2, 0], [0, 0]]],
     ((-3 + cos(theta)**2 - cos(2 * phi) * sin(theta)**2) * -3 / 8)],
    ["0+ → 1+ → 2+ → 0+", [theta, phi, [0, 1], [1, 1, 0], [[2, 0], [0, 0]]],
     ((-3 + cos(theta)**2 + cos(2 * phi) * sin(theta)**2) * -3 / 8)],
    ["0+ → 1- → 2+ → 2+", [theta, phi, [0, 1], [1, 0, 0], [[2, 0], [2, 0]]],
     ((-29 + 7 * cos(theta)**2 - 7 * cos(2 * phi) * sin(theta)**2) * -3 / 80)],
    ])
def test_double(description, W_args, result):
    print(f"Calculating angular distribution {description}.")
    assert is_equal(W(*W_args), result)


@pytest.mark.parametrize("description, W_args, result", [
    ["0+ → 1- → 0+ → 2+ → 0+",
     [theta, phi, [0, 1], [1, 0, 0], [[0, 0], [2, 0], [0, 0]]], 1],
    ["0+ → 1- → 2+ → 2+ → 0+",
     [theta, phi, [0, 1], [1, 0, 0], [[2, 0], [2, 0], [0, 0]]],
     ((17 - 3 * cos(theta)**2 + 3 * cos(2 * phi) * sin(theta)**2) / 16)],
    ["0+ → 1+ → 2+ → 2+ → 0+",
     [theta, phi, [0, 1], [1, 1, 0], [[2, 0], [2, 0], [0, 0]]],
     ((17 - 3 * cos(theta)**2 - 3 * cos(2 * phi) * sin(theta)**2) / 16)],
    ])
def test_triple(description, W_args, result):
    print(f"Calculating angular distribution {description}.")
    assert is_equal(W(*W_args), result)


@pytest.mark.parametrize("description, W_args, result", [
    ["0+ → 1- → 0+",
     [theta, phi, [0, 1], [1, 0, 0], [[0, delta]]],
     ((3 * (1 + cos(theta)**2 - cos(2 * phi) *
         sin(theta)**2))/(4 * (1 + delta**2)))],
     ])
def test_delta(description, W_args, result):
    print(f"Calculating angular distribution {description}.")
    assert is_equal(W(*W_args), result)
