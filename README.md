# Entropy-Stable Discontinuous Spectral-Element Methods for the Spherical Shallow Water Equations in Covariant Form

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT) [![doi: https://doi.org/10.5281/zenodo.17059901](https://zenodo.org/badge/DOI/10.5281/zenodo.17059901.svg)](https://doi.org/10.5281/zenodo.17059901)

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

Please cite the above manuscript if you use or adapt the code in this repository or that of the associated solvers implemented in [TrixiAtmo.jl](https://github.com/trixi-framework/TrixiAtmo.jl) in your research. Any questions regarding the content of the manuscript or technical issues with the code should be directed to the corresponding author at montoya.tristan@gmail.com.

## Abstract
We introduce discontinuous spectral-element methods of arbitrary order that are well balanced, conservative of mass, and conservative or dissipative of total energy (i.e., a mathematical entropy function) for the rotating shallow water equations on curved manifolds such as the sphere. The proposed methods are based on a skew-symmetric flux formulation in which the prognostic variables correspond to the fluid layer depth and contravariant momentum components, which we implement and analyze within a general flux-differencing framework using tensor-product summation-by-parts operators. Such schemes are shown to guarantee mass and energy conservation (or energy dissipation for an appropriate choice of numerical flux at element interfaces) at the semi-discrete level on general unstructured quadrilateral grids, in addition to well balancing for arbitrary continuous bottom topographies. Furthermore, our semi-discrete covariant formulation permits an analytical representation of the geometry and associated metric terms while satisfying the aforementioned entropy stability, conservation, and well-balancing properties without the need to approximate the metric terms so as to satisfy discrete metric identities. Numerical experiments are presented in order to verify the above structure preservation properties as well as to assess the accuracy and robustness of the proposed schemes within the context of several standard test cases characteristic of idealized atmospheric flows.

## Installation
First, make you have [Julia](https://julialang.org/downloads/) installed you haven't already done so. Then, assuming that you are using Linux or macOS and have git installed, follow the steps below.

1. Clone this repository by entering the command `git clone https://github.com/tristanmontoya/paper-2025-spherical-shallow-water.git` in the terminal.

2. Within the `code` directory, use the command `julia --project=.` to open the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) and activate the project within the current directory. 

3. Install all dependencies by entering `using Pkg; Pkg.instantiate()` in the REPL.

## Reproducibility instructions
Here, we describe how to generate the results using the provided drivers, and how to produce the results in the manuscript using the provided Jupyter notebooks. Note that the tests run significantly faster with multithreading enabled (for example, add `--threads 8` to the `julia` command if you want to use eight threads) as this allows for local element-based operations to be performed simultaneously. If using multiple Julia threads, it is [usually best to set the number of BLAS threads to 1](https://carstenbauer.github.io/ThreadPinning.jl/stable/examples/ex_blas/#Beware:-Interaction-between-Julia-threads-and-BLAS-threads) (for example, using the `OPENBLAS_NUM_THREADS` environment variable).

To run the four test cases associated with Sections 5.1, 5.2, 5.3, and 5.4, the following commands should be run in the Julia REPL after following the installation instructions.

```julia
using Trixi, TrixiAtmo, SphericalShallowWater
run_unsteady_solid_body_rotation()  # Section 5.1 
run_isolated_mountain()  # Section 5.2
run_barotropic_instability()  # Section 5.3
run_rossby_haurwitz()  # Section 5.4
```

This will automatically generate the `results/` directory, where the subdirectories for each run will be organized by date, test case name, and an additional descriptive suffix for example, `20250713_isolated_mountain_ec/`. Depending on the computational resources you have available, it may be helpful to break up the calls to `run_driver` within the above function calls in order to run more cases in parallel.

Reproducing the plots in Figures 2, 3, 5, and 8 is similarly automated through the following REPL commands, which will place `.pdf` files generated with [Makie.jl](https://github.com/MakieOrg/Makie.jl/) in the (automatically generated) `plots/` directory:

```julia
using Trixi, TrixiAtmo, SphericalShallowWater
plot_unsteady_solid_body_rotation()  # Figure 2
plot_isolated_mountain()  # Figure 3
plot_barotropic_instability()  # Figure 5
plot_rossby_haurwitz()  # Figure 8
```

 Figures 4, 7, and 8 are contour plots of relative vorticity for the isolated mountain case as well as the barotropic instability with and without perturbation. After using [Trixi2Vtk.jl](https://github.com/trixi-framework/Trixi2Vtk.jl) to convert the `.h5` files generated by Trixi.jl to `.vtu` files, the provided state files `paraview_states/isolated_mountain.pvsm` (Figure 4) and `paraview_states/barotropic_instability.pvsm` (Figures 7 and 8) can used to visualize the results in [ParaView](https://www.paraview.org/).

## License
This repository is released under the [MIT license](https://github.com/tristanmontoya/paper-2025-spherical-shallow-water/blob/main/LICENSE).

