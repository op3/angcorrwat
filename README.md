# angcorrwat

This python library can be used to calculate angular distributions.
It is planned to add the calculation of angular correlations in the future.


## Usage

To calculate the ϑ-, φ- and δ-dependent angular distribution of a 0⁺ → 1⁻ → 0 cascade, use the following code:

```
    from sympy import Symbol
    from angcorrwat.angdist import W

    theta = Symbol('theta')
    phi = Symbol('phi')
    delta = Symbol('delta')

    W(theta, phi, [0, 1], [[1, 0, 0], [0, 1, delta]])
```

The incident γ-ray beam is assumed to be linearily polarized.
The resulting angular distribution is normalized to 1.


## License<a name="license"></a>

Copyright © 2018

O. Papst (opapst@ikp.tu-darmstadt.de)

angcorrwat is distributed under the terms of the GNU General Public License, version 3 or later.
See the `LICENSE` file.


## References

<a name="ref-1">[1]</a> R. M. Steffen *et al.*, “Angular distribution and correlation of gamma rays”, in *The electromagnetic interaction in nuclear spectroscopy*, edited by W. D. Hamilton (North-Holland, Amsterdam, 1975) Chap. 12, pp. 505–582, ISBN: 978-0-4441-0519-6.

<a name="ref-2">[2]</a> K. S. *Krane et al.*, “Directional correlations of gamma radiations emitted from nuclear states oriented by nuclear reactions or cryogenic methods”, At. Data Nucl. Data Tables **11**, 351 (1973). [`doi:10.1016/S0092-640X(73)80016-6`](https://doi.org/10.1016/S0092-640X(73)80016-6).  
