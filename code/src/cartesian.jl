# Initial condition is atmosphere at rest with constant total geopotential height
function initial_condition_well_balanced(x, t, equations::ShallowWaterEquations3D)
    # Constant total geopotential height
    H = 5960.0
    v1 = v2 = v3 = 0.0

    # Non-constant topography
    b = bottom_topography_isolated_mountain(x)

    return SVector(H, v1, v2, v3, b)
end

function Trixi.calc_error_norms(
    func,
    u,
    t,
    analyzer,
    mesh::Union{
        StructuredMesh{2},
        StructuredMeshView{2},
        UnstructuredMesh2D,
        P4estMesh{2},
        P4estMeshView{2},
        T8codeMesh{2},
    },
    equations::ShallowWaterEquations3D,
    initial_condition,
    dg::DGSEM,
    cache,
    cache_analysis,
)
    (; vandermonde, weights) = analyzer
    (; node_coordinates, inverse_jacobian) = cache.elements
    (; u_local, u_tmp1, x_local, x_tmp1, jacobian_local, jacobian_tmp1) = cache_analysis

    # Set up data structures
    l2_error = zero(func(Trixi.get_node_vars(u, equations, dg, 1, 1, 1), equations))
    linf_error = copy(l2_error)
    l2_normalization = copy(l2_error)
    linf_normalization = copy(l2_error)

    # Iterate over all elements for error calculations
    for element in Trixi.eachelement(dg, cache)
        # Interpolate solution and node locations to analysis nodes
        Trixi.multiply_dimensionwise!(
            u_local,
            vandermonde,
            view(u, :, :, :, element),
            u_tmp1,
        )
        Trixi.multiply_dimensionwise!(
            x_local,
            vandermonde,
            view(node_coordinates, :, :, :, element),
            x_tmp1,
        )
        Trixi.multiply_scalar_dimensionwise!(
            jacobian_local,
            vandermonde,
            inv.(view(inverse_jacobian, :, :, element)),
            jacobian_tmp1,
        )

        # Calculate errors at each analysis node
        for j in Trixi.eachnode(analyzer), i in Trixi.eachnode(analyzer)
            u_exact = initial_condition(
                Trixi.get_node_coords(x_local, equations, dg, i, j),
                t,
                equations,
            )
            func_exact = func(Trixi.get_node_vars(u_local, equations, dg, i, j), equations)

            diff = func(u_exact, equations) - func_exact
            # We take absolute value as we need the Jacobian here for the volume calculation
            abs_jacobian_local_ij = abs(jacobian_local[i, j])

            l2_error += diff .^ 2 * (weights[i] * weights[j] * abs_jacobian_local_ij)
            linf_error = @. max(linf_error, abs(diff))

            linf_normalization = @. max(linf_normalization, abs(func_exact))
            l2_normalization +=
                func_exact .^ 2 * (weights[i] * weights[j] * abs_jacobian_local_ij)
        end
    end

    # Normalize both errors 
    l2_error = @. sqrt(l2_error / l2_normalization)
    linf_error = @. linf_error / linf_normalization

    return l2_error, linf_error
end

function run_isolated_mountain_cartesian()
    run_driver(
        "elixirs/elixir_spherical_shallow_water_cartesian.jl",
        1,
        polydeg = 3,
        initial_condition = initial_condition_well_balanced,
        surface_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal),
        initial_cells_per_dimension = 20,
        identifier = "_cartesian_ec",
        interval = 50,
        tspan = (0.0, 15.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 150,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water_cartesian.jl",
        1,
        polydeg = 3,
        initial_condition = initial_condition_well_balanced,
        surface_flux = (
            FluxPlusDissipation(
                flux_wintermeyer_etal,
                DissipationLaxFriedrichsEntropyVariables(),
            ),
            flux_nonconservative_wintermeyer_etal,
        ),
        initial_cells_per_dimension = 20,
        identifier = "_cartesian_es",
        interval = 50,
        tspan = (0.0, 15.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 150,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water_cartesian.jl",
        1,
        polydeg = 3,
        initial_condition = initial_condition_isolated_mountain,
        surface_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal),
        initial_cells_per_dimension = 20,
        identifier = "_cartesian_ec",
        interval = 50,
        tspan = (0.0, 15.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 150,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water_cartesian.jl",
        1,
        polydeg = 3,
        initial_condition = initial_condition_isolated_mountain,
        surface_flux = (
            FluxPlusDissipation(
                flux_wintermeyer_etal,
                DissipationLaxFriedrichsEntropyVariables(),
            ),
            flux_nonconservative_wintermeyer_etal,
        ),
        initial_cells_per_dimension = 20,
        identifier = "_cartesian_es",
        interval = 50,
        tspan = (0.0, 15.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 150,
    )
end

function plot_unsteady_solid_body_rotation_cartesian()
    # Figure 2a
    plot_convergence(
        [
            "../results/20250714_unsteady_solid_body_rotation_N3_cartesian_ec/",
            "../results/20250714_unsteady_solid_body_rotation_N3_cartesian_es/",
        ],
        plot_name = "unsteady_solid_body_rotation_convergence_N3_cartesian.pdf",
        triangle_bottom_order = 4,
        triangle_top_order = 3,
    )
    # Figure 2b
    plot_convergence(
        [
            "../results/20250714_unsteady_solid_body_rotation_N4_cartesian_ec/",
            "../results/20250714_unsteady_solid_body_rotation_N4_cartesian_es/",
        ],
        plot_name = "unsteady_solid_body_rotation_convergence_N4_cartesian.pdf",
        triangle_bottom_order = 5,
        triangle_top_order = 5,
    )

    # Figure 2c
    plot_convergence(
        [
            "../results/20250714_unsteady_solid_body_rotation_p_refine_cartesian_ec/",
            "../results/20250714_unsteady_solid_body_rotation_p_refine_cartesian_es/",
        ],
        triangle_bottom = false,
        triangle_top = false,
        plot_name = "unsteady_solid_body_rotation_p_refine_M4_cartesian.pdf",
        xkey = "N",
        xlabel = L"Polynomial degree $N$",
        xticks = collect(2:10),
        xscale = identity,
        legend_position = (:left, :bottom),
    )
end

function plot_isolated_mountain_cartesian()
    # Figure 3a
    plot_evolution(
        [
            "../results/20250715_well_balanced_cartesian_ec/N3M20",
            "../results/20250715_well_balanced_cartesian_es/N3M20",
        ],
        plot_name = "well_balanced_l2_h_evolution_N3M20_cartesian.pdf",
        ykey = "l2_h",
        legend_position = (:left, :top),
        relative = false,
        ylabel = L"Normalized $L^2$ height error",
        xlims = [0, 15],
        ylims= [-0.2,5],
        xticks = [0, 5, 10, 15],
        ynorm = 1e-14,
        exponent_text = L"\times 10^{-14}",
    )

    # Figure 3b
    plot_evolution(
        [
            "../results/20250715_isolated_mountain_cartesian_ec/N3M20",
            "../results/20250715_isolated_mountain_cartesian_es/N3M20",
        ],
        plot_name = "isolated_mountain_mass_evolution_N3M20_cartesian.pdf",
        legend_position = (:left, :top),
        xlims = [0, 15],
        xticks = [0, 5, 10, 15],
        ylims= [-0.85,4.5],
        ykey = "waterheight",
        ylabel = LaTeXString("Normalized mass change"),
        ynorm = 1e-14,
        exponent_text = L"\times 10^{-14}",
    )

    # Figure 3c
    plot_evolution(
        [
            "../results/20250715_isolated_mountain_cartesian_ec/N3M20",
            "../results/20250715_isolated_mountain_cartesian_es/N3M20",
        ],
        plot_name = "isolated_mountain_entropy_evolution_N3M20_cartesian.pdf",
        legend_position = (:left, :bottom),
        xlims = [0, 15],
        xticks = [0, 5, 10, 15],
        ylims = [-8.25, 0.25],
        ykey = "e_total",
        ynorm = 1e-8,
        exponent_text = L"\times 10^{-8}",
    )
end
