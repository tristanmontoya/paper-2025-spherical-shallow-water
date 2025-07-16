#!/bin/bash

nohup julia --threads 20 --project=. -e 'using SphericalShallowWater, TrixiAtmo, Trixi;  run_driver(
        "elixirs/elixir_spherical_shallow_water_cartesian.jl",
        5,
        polydeg = 3,
        initial_condition = initial_condition_unsteady_solid_body_rotation,
        surface_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal),
        initial_cells_per_dimension = 4,
        identifier = "_N3_cartesian_ec",
        interval = 50,
        tspan = (0.0, 5.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 50,
    )' > p3_ec.txt &

nohup julia --threads 20 --project=. -e 'using SphericalShallowWater, TrixiAtmo, Trixi;  run_driver(
        "elixirs/elixir_spherical_shallow_water_cartesian.jl",
        5,
        polydeg = 3,
        initial_condition = initial_condition_unsteady_solid_body_rotation,
        surface_flux = (
            FluxPlusDissipation(
                flux_wintermeyer_etal,
                DissipationLaxFriedrichsEntropyVariables(),
            ),
            flux_nonconservative_wintermeyer_etal,
        ),
        initial_cells_per_dimension = 4,
        identifier = "_N3_cartesian_es",
        interval = 50,
	tspan = (0.0, 5.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 50,
    )' > p3_es.txt &

nohup julia --threads 20 --project=. -e 'using SphericalShallowWater, TrixiAtmo, Trixi;  run_driver(
        "elixirs/elixir_spherical_shallow_water_cartesian.jl",
        5,
        polydeg = 4,
        initial_condition = initial_condition_unsteady_solid_body_rotation,
        surface_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal),
        initial_cells_per_dimension = 4,
        identifier = "_N4_cartesian_ec",
        interval = 50,
        tspan = (0.0, 5.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 50,
    )' > p4_ec.txt &


nohup julia --threads 20 --project=. -e 'using SphericalShallowWater, TrixiAtmo, Trixi;  run_driver(
        "elixirs/elixir_spherical_shallow_water_cartesian.jl",
        5,
        polydeg = 4,
        initial_condition = initial_condition_unsteady_solid_body_rotation,
        surface_flux = (
            FluxPlusDissipation(
                flux_wintermeyer_etal,
                DissipationLaxFriedrichsEntropyVariables(),
            ),
            flux_nonconservative_wintermeyer_etal,
        ),
        initial_cells_per_dimension = 4,
        identifier = "_N4_cartesian_es",
        interval = 50,
	tspan = (0.0, 5.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 50,
    )' > p4_es.txt &
