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

## Installation
First, make you have [Julia](https://julialang.org/downloads/) installed you haven't already done so (we recommend using the latest stable release). Then, assuming that you are using Linux or macOS and have git installed, follow the steps below.

1. Clone this repository by entering the command `git clone https://github.com/tristanmontoya/paper-2025-spherical-shallow-water.git` in the terminal.

2. Within the `code` directory, use the command `julia --project=.` to open the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) and activate the project within the current directory. 

3. Install all dependencies by entering `using Pkg; Pkg.instantiate()` in the REPL.

## Reproducibility instructions
Here, we describe how to generate the results using the provided drivers, and how to produce the results in the manuscript using the provided Jupyter notebooks. Note that the tests run significantly faster with multithreading enabled (for example, add `--threads 8` to the `julia` command if you want to use eight threads) as this allows for local element-based operations to be performed simultaneously. If using multiple Julia threads, it is [usually best to set the number of BLAS threads to 1](https://carstenbauer.github.io/ThreadPinning.jl/dev/explanations/blas/) (for example, using the `OPENBLAS_NUM_THREADS` environment variable). For all tests, enter the top-level directory and run `julia --project=.` and then load the driver package by entering `using SphericalShallowWater` in the REPL.
