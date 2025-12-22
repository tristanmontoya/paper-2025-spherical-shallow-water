function run_unsteady_solid_body_rotation()
    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        5,
        polydeg = 3,
        initial_condition = initial_condition_unsteady_solid_body_rotation,
        auxiliary_field = bottom_topography_unsteady_solid_body_rotation,
        surface_flux = surface_flux_ec,
        initial_cells_per_dimension = 4,
        identifier = "_ec_N3",
        interval = 50,
        tspan = (0.0, 5.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 50,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        5,
        polydeg = 3,
        initial_condition = initial_condition_unsteady_solid_body_rotation,
        auxiliary_field = bottom_topography_unsteady_solid_body_rotation,
        surface_flux = surface_flux_es,
        initial_cells_per_dimension = 4,
        identifier = "_es_N3",
        interval = 50,
        tspan = (0.0, 5.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 50,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        5,
        polydeg = 4,
        initial_condition = initial_condition_unsteady_solid_body_rotation,
        auxiliary_field = bottom_topography_unsteady_solid_body_rotation,
        surface_flux = surface_flux_ec,
        initial_cells_per_dimension = 4,
        identifier = "_ec_N4",
        interval = 50,
        tspan = (0.0, 5.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 50,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        5,
        polydeg = 4,
        initial_condition = initial_condition_unsteady_solid_body_rotation,
        auxiliary_field = bottom_topography_unsteady_solid_body_rotation,
        surface_flux = surface_flux_es,
        initial_cells_per_dimension = 4,
        identifier = "_es_N4",
        interval = 50,
        tspan = (0.0, 5.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 50,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        1,
        polydeg = 2:10,
        initial_condition = initial_condition_unsteady_solid_body_rotation,
        auxiliary_field = bottom_topography_unsteady_solid_body_rotation,
        surface_flux = surface_flux_ec,
        initial_cells_per_dimension = 4,
        identifier = "_ec_p_refine",
        interval = 50,
        tspan = (0.0, 5.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 50,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        1,
        polydeg = 2:10,
        initial_condition = initial_condition_unsteady_solid_body_rotation,
        auxiliary_field = bottom_topography_unsteady_solid_body_rotation,
        surface_flux = surface_flux_es,
        initial_cells_per_dimension = 4,
        identifier = "_es_p_refine",
        interval = 50,
        tspan = (0.0, 5.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 50,
    )
end

function run_isolated_mountain()
    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        1,
        polydeg = 3,
        initial_condition = initial_condition_well_balanced,
        auxiliary_field = bottom_topography_isolated_mountain,
        surface_flux = surface_flux_ec,
        initial_cells_per_dimension = 20,
        identifier = "_ec",
        interval = 50,
        tspan = (0.0, 15.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 150,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        1,
        polydeg = 3,
        initial_condition = initial_condition_well_balanced,
        auxiliary_field = bottom_topography_isolated_mountain,
        surface_flux = surface_flux_es,
        initial_cells_per_dimension = 20,
        identifier = "_es",
        interval = 50,
        tspan = (0.0, 15.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 150,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        1,
        polydeg = 3,
        initial_condition = initial_condition_isolated_mountain,
        auxiliary_field = bottom_topography_isolated_mountain,
        surface_flux = surface_flux_ec,
        initial_cells_per_dimension = 20,
        identifier = "_ec",
        interval = 50,
        tspan = (0.0, 15.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 150,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        1,
        polydeg = 3,
        initial_condition = initial_condition_isolated_mountain,
        auxiliary_field = bottom_topography_isolated_mountain,
        surface_flux = surface_flux_es,
        initial_cells_per_dimension = 20,
        identifier = "_es",
        interval = 50,
        tspan = (0.0, 15.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 150,
    )
end

function run_barotropic_instability()
    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        3,
        polydeg = 3,
        initial_condition = initial_condition_steady_barotropic_instability,
        auxiliary_field = nothing,
        surface_flux = surface_flux_ec,
        initial_cells_per_dimension = 16,
        identifier = "_ec",
        interval = 50,
        tspan = (0.0, 12.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 120,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        3,
        polydeg = 3,
        initial_condition = initial_condition_steady_barotropic_instability,
        auxiliary_field = nothing,
        surface_flux = surface_flux_es,
        initial_cells_per_dimension = 16,
        identifier = "_es",
        interval = 50,
        tspan = (0.0, 12.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 120,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        3,
        polydeg = 3,
        initial_condition = initial_condition_barotropic_instability,
        auxiliary_field = nothing,
        surface_flux = surface_flux_ec,
        initial_cells_per_dimension = 16,
        identifier = "_ec",
        interval = 50,
        tspan = (0.0, 12.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 120,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        3,
        polydeg = 3,
        initial_condition = initial_condition_barotropic_instability,
        auxiliary_field = nothing,
        surface_flux = surface_flux_es,
        initial_cells_per_dimension = 16,
        identifier = "_es",
        interval = 50,
        tspan = (0.0, 12.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 120,
    )

    run_timestep_study( # New for revised version
        "elixirs/elixir_spherical_shallow_water.jl",
        6,
        Float64;
        results_dir = RESULTS_DIR,
        polydeg = 3,
        cells_per_dimension = 16,
        n_saves = 120,
        initial_condition = initial_condition_steady_barotropic_instability,
        auxiliary_field = nothing,
        surface_flux = surface_flux_ec,
        tspan = (0.0, 12.0 * SECONDS_PER_DAY),
        initial_dt = 320.0,
    )

    run_driver( # New for revised version
           "elixirs/elixir_spherical_shallow_water_standard_dg.jl",
           1,
           polydeg = 3,
           initial_condition = initial_condition_barotropic_instability,
           auxiliary_field = nothing,
           initial_cells_per_dimension = 64,
           identifier = "_standard",
           interval = 50,
           tspan = (0.0, 12.0 * SECONDS_PER_DAY),
           cfl = 0.1,
           n_saves = 120,
       )
end

function run_rossby_haurwitz()
    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        1,
        polydeg = 3,
        initial_condition = initial_condition_rossby_haurwitz,
        auxiliary_field = nothing,
        surface_flux = surface_flux_ec,
        initial_cells_per_dimension = 16,
        identifier = "_ec_N3",
        interval = 50,
        tspan = (0.0, 28.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 280,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        1,
        polydeg = 3,
        initial_condition = initial_condition_rossby_haurwitz,
        auxiliary_field = nothing,
        surface_flux = surface_flux_es,
        initial_cells_per_dimension = 16,
        identifier = "_es_N3",
        interval = 50,
        tspan = (0.0, 28.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 280,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water_standard_dg.jl",
        1,
        polydeg = 3,
        initial_condition = initial_condition_rossby_haurwitz,
        auxiliary_field = nothing,
        initial_cells_per_dimension = 16,
        identifier = "_standard_N3",
        interval = 50,
        tspan = (0.0, 28.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 280,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        1,
        polydeg = 6,
        initial_condition = initial_condition_rossby_haurwitz,
        auxiliary_field = nothing,
        surface_flux = surface_flux_ec,
        initial_cells_per_dimension = 8,
        identifier = "_ec_N6",
        interval = 50,
        tspan = (0.0, 28.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 280,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        1,
        polydeg = 6,
        initial_condition = initial_condition_rossby_haurwitz,
        auxiliary_field = nothing,
        surface_flux = surface_flux_es,
        initial_cells_per_dimension = 8,
        identifier = "_es_N6",
        interval = 50,
        tspan = (0.0, 28.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 280,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water_standard_dg.jl",
        1,
        polydeg = 6,
        initial_condition = initial_condition_rossby_haurwitz,
        auxiliary_field = nothing,
        initial_cells_per_dimension = 8,
        identifier = "_standard_N6",
        interval = 50,
        tspan = (0.0, 28.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 280,
    )
end