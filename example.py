#!/usr/bin/env python

from sympy import Symbol
from angdist import W

theta = Symbol('theta')
phi = Symbol('phi')

print("Angular distribution of 0⁺ → 1⁻ → 0⁺")
#res = W(theta, phi, 0, 1, 0, 0, 1, 1, 2, 0, 1, 2, 0)
res = W(theta, phi, 0, 0, 1, 1, 2, 0, [[1, 0, 1, 2]], 0)
print(res)

print("Angular distribution of 0⁺ → 1⁺ → 0⁺")
#res = W(theta, phi, 0, 1, 0, 1, 1, 0, 2, 0, 1, 2, 0)
res = W(theta, phi, 0, 1, 1, 0, 2, 0, [[1, 0, 1, 2]], 0)
print(res)

print("Angular distribution of 0⁺ → 1⁺ → 2⁺")
#res = W(theta, phi, 0, 1, 2, 1, 1, 0, 2, 0, 1, 2, 0)
res = W(theta, phi, 0, 1, 1, 0, 2, 0, [[1, 1, 2, 0]], 2)
print(res)

print("Angular distribution of 0⁺ → 1⁻ → 2⁺")
#res = W(theta, phi, 0, 1, 2, 0, 1, 1, 2, 0, 1, 2, 0)
res = W(theta, phi, 0, 0, 1, 1, 2, 0, [[1, 1, 2, 0]], 1)
print(res)

print("Angular distribution of 0⁺ → 2⁺ → 0⁺")
#res = W(theta, phi, 0, 2, 0, 0, 2, 1, 3, 0, 2, 3, 0)
res = W(theta, phi, 0, 0, 2, 1, 3, 0, [[2, 2, 3, 0]], 0)
print(res)

print("Wrong:")
print("Angular distribution of 0⁺ → 1⁻ → 2⁺ → 0⁺")
res = W(theta, phi, 0, 0, 1, 1, 2, 0, [[1, 1, 2, 0], [2, 2, 3, 0]], 2)
print(res)
