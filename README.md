# angcorrwat

[![Build Status](https://travis-ci.com/op3/angcorrwat.svg?branch=master)](https://travis-ci.com/op3/angcorrwat)
[![License](https://img.shields.io/badge/License-GPL%20v3+-blue.svg)](LICENSE)
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/op3/angcorrwat/master?filepath=doc%2Ftutorial.ipynb)

This python library can be used to calculate angular distributions/correlations.


## Dependencies

To use the library, the following dependencies have to be fullfilled:

- python≥3.6
- sympy

To run the tutorial notebook, matplotlib is required.

## Usage

To calculate the ϑ-, φ- and δ-dependent angular distribution of a 0⁺ → 1⁻ → 0 cascade, use the following code:

```
from sympy import Symbol
from angcorrwat import W

theta = Symbol('theta')
phi = Symbol('phi')
delta = Symbol('delta')

W(theta, phi, [0, 1], [1, 0, 0], [[0, delta]])
```

The resulting angular distribution is normalized to 1, when integrated over the solid angle.
The incident γ-ray beam is assumed to be linearily polarized.
Spherical coordinates are used:
While the azimuthal angle φ refers to the angle between the polarization plane and the reaction plane,
the polar angle ϑ is given with respect to the incident γ-ray.
For the multipole mixing ratio δ, the convention by Krane, Steffen and Wheeler [\[2\]](#ref-2) is used.

The arguments `theta` and `phi` refer to the respective spherical coordinates ϑ and φ.
The third argument `[0, 1]` = `[J, π]` refers to the initial state J<sup>π</sup> = 0<sup>+</sup> (0 is negative and 1 positive parity).
The fourth argument `[1, 0, 0]` = `[J, π, δ]` refers to the state (de)excited by the first γ, with `δ` refering to the multipole mixing ratio for the transition between the initial and excited state.
Finally, a list of states `[[J, δ], …]` of the subsequent cascade is given, with `δ` refering to the multipole mixing ratio for the transition between the previous and current state.
The angular distribution of the final state in the cascade is returned.

For more examples, take a look at the [tutorial](doc/tutorial.ipynb) ([launch in binder](https://mybinder.org/v2/gh/op3/angcorrwat/master?filepath=doc%2Ftutorial.ipynb)).

## License<a name="license"></a>

© 2018 O. Papst [`<opapst@ikp.tu-darmstadt.de>`](mailto:opapst@ikp.tu-darmstadt.de)

angcorrwat is distributed under the terms of the GNU General Public License, version 3 or later.
See the [`LICENSE`](LICENSE) file.


## References

<a name="ref-1">[1]</a> R. M. Steffen *et al.*, “Angular distribution and correlation of gamma rays”, in *The electromagnetic interaction in nuclear spectroscopy*, edited by W. D. Hamilton (North-Holland, Amsterdam, 1975) Chap. 12, pp. 505–582, ISBN: 978-0-4441-0519-6.

<a name="ref-2">[2]</a> K. S. Krane *et al.*, “Directional correlations of gamma radiations emitted from nuclear states oriented by nuclear reactions or cryogenic methods”, At. Data Nucl. Data Tables **11**, 351 (1973). [`doi:10.1016/S0092-640X(73)80016-6`](https://doi.org/10.1016/S0092-640X(73)80016-6).  
