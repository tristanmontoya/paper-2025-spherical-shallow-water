# Entropy-Stable Discontinuous Spectral-Element Methods for the Spherical Shallow Water Equations in Covariant Form

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

This repository contains information and code to reproduce the results presented in the following article:
```bibtex
@online{MontoyaSphericalShallowWater2025,
  title={Entropy-Stable Discontinuous Spectral-Element Methods for the 
         Spherical Shallow Water Equations in Covariant Form},
  author={Montoya, Tristan and Rueda-Ramírez, Andrés M. and Gassner, Gregor J.},
  year={2025},
  journal={In preparation}
}
```

## Abstract
We introduce discontinuous spectral-element methods of arbitrary order for general unstructured quadrilateral grids which are well balanced, conservative of mass, and conservative or dissipative of total energy (i.e. a mathematical entropy function) for the rotating shallow water equations on curved manifolds such as the sphere. The proposed methods are based on a skew-symmetric flux formulation in which the prognostic variables correspond to the geopotential height and contravariant momentum components, which we implement and analyze within a general flux-differencing framework using tensor-product summation-by-parts operators. Such schemes are shown to guarantee mass and energy conservation (or energy dissipation for an appropriate choice of numerical flux at element interfaces) in addition to well-balancedness for arbitrary continuous bottom topographies, and do not require any discrete metric identities to be satisfied, thereby permitting an analytical representation of the geometry and associated metric terms. Numerical experiments are presented in order to verify the above structure-preserving properties as well as to assess the accuracy and robustness of the proposed schemes within the context of several standard test cases characteristic of idealized atmospheric flows.