module SphericalShallowWater

using Trixi, TrixiAtmo, Trixi2Vtk
using CairoMakie, LaTeXStrings, Dates, Printf, CSV

export EXAMPLES_DIR, RESULTS_DIR
export surface_flux_ec, surface_flux_es
export run_driver, plot_convergence, plot_evolution, calc_norms
export run_unsteady_solid_body_rotation,
    run_isolated_mountain, run_barotropic_instability, run_rossby_haurwitz
export plot_unsteady_solid_body_rotation,
    plot_isolated_mountain, plot_barotropic_instability, plot_rossby_haurwitz
export initial_condition_well_balanced, initial_condition_steady_barotropic_instability

const EXAMPLES_DIR = TrixiAtmo.examples_dir()
const RESULTS_DIR = joinpath(dirname(dirname(@__DIR__)), "results")
const PLOTS_DIR = joinpath(dirname(dirname(@__DIR__)), "plots")

# EC and ES surface fluxes
const surface_flux_ec = (flux_ec, flux_nonconservative_surface_simplified)
const surface_flux_es = (
    FluxPlusDissipation(flux_ec, DissipationLocalLaxFriedrichs()),
    flux_nonconservative_surface_simplified,
)

function run_driver(
    elixir::AbstractString,
    iterations = 1, # number of times to double `cells_per_dimension`
    RealT = Float64;
    results_dir = RESULTS_DIR, # top level directory to store results
    polydeg = 3, # can also be a single integer value
    initial_cells_per_dimension = 4, # initial value of `cells_per_dimension`
    initial_condition = initial_condition_unsteady_solid_body_rotation,
    date = Dates.format(today(), dateformat"yyyymmdd"),
    identifier = "", # suffix to identify a specific run
    kwargs...,
)

    mod = @__MODULE__

    # Get name for project
    ic_name = replace(string(initial_condition), r"^initial_condition_" => "")
    project_dir = mkpath(joinpath(results_dir, string(date, "_", ic_name, identifier)))
    println("Project directory: ", project_dir)

    # Format top-level output file and write headers
    fmt1 = Printf.Format(
        "%-4d" * "%-4d" * "%-2.5f " * "%-25.17e"^5 * "missing  " * "missing" * "\n",
    ) # no EOC
    fmt2 = Printf.Format(
        "%-4d" * "%-4d" * "%-2.5f " * "%-25.17e"^5 * "%-9.2f" * "%-9.2f" * "\n",
    )
    headers = [
        "N   ",
        "M   ",
        "end_time ",
        "resolution_km",
        "l2_depth_error",
        "linf_depth_error",
        "l2_height_error",
        "linf_height_error",
        "l2_order ",
        "linf_order",
    ]
    open(joinpath(project_dir, "analysis.dat"), "w") do io
        println(
            io,
            string(
                headers[1:3]...,
                rpad.(headers[4:end-2], 25)...,
                headers[end-1],
                headers[end],
            ),
        )
    end

    resolutions = RealT[]
    l2_errors = RealT[]
    linf_errors = RealT[]
    end_times = RealT[]

    cells_per_dimension = initial_cells_per_dimension .* 2 .^ ((1:iterations) .- 1)

    # run simulations and extract errors
    for N in polydeg
        for M in cells_per_dimension
            trixi_include(
                mod,
                elixir;
                kwargs...,
                output_dir = joinpath(project_dir, string("N", N, "M", M)),
                initial_condition = initial_condition,
                polydeg = N,
                cells_per_dimension = M,
            )

            resolution_km = π * EARTH_RADIUS / (2 * M * N * 1000) # scaled by 1000 to get km
            end_time = mod.t_final / SECONDS_PER_DAY

            append!(end_times, end_time)
            append!(resolutions, resolution_km)
            append!(l2_errors, mod.l2_height_error)
            append!(linf_errors, mod.linf_height_error)

            open(joinpath(project_dir, "analysis.dat"), "a") do io
                if M == initial_cells_per_dimension # don't compute order for first grid
                    Printf.format(
                        io,
                        fmt1,
                        N,
                        M,
                        end_time,
                        resolution_km,
                        mod.l2_depth_error,
                        mod.linf_depth_error,
                        l2_height_error,
                        linf_height_error,
                    )
                else
                    Printf.format(
                        io,
                        fmt2,
                        N,
                        M,
                        end_time,
                        resolution_km,
                        mod.l2_depth_error,
                        mod.linf_depth_error,
                        mod.l2_height_error,
                        mod.linf_height_error,
                        log(l2_errors[end] / l2_errors[end-1]) / log(1 / 2),
                        log(linf_errors[end] / linf_errors[end-1]) / log(1 / 2),
                    )
                end
            end

            println("\n\n")
            println("#"^100)
        end
    end

end

function plot_convergence(
    dirs = joinpath(
        RESULTS_DIR,
        string(
            Dates.format(today(), dateformat"yyyymmdd"),
            "_unsteady_solid_body_rotation",
        ),
    );
    plots_dir = PLOTS_DIR,
    plot_name = nothing,
    labels = [LaTeXString("EC"), LaTeXString("ES")],
    styles = [:dash, :solid],
    colors = 1:length(labels),
    file = "analysis.dat",
    xkey = "resolution_km",
    ykey = "l2_height_error",
    ynorm = nothing,
    legend = true,
    xlabel = LaTeXString("Nominal resolution (km)"),
    ylabel = L"Normalized $L^2$ height error",
    font = "CMU Serif",
    size = (350, 350),
    fontsize = 15,
    legendfontsize = 15,
    trianglefontsize = 13,
    linewidth = 1.5,
    xticklabelfont = "CMU Serif",
    yticklabelfont = "CMU Serif",
    legendfont = "CMU Serif",
    xlims = nothing,
    ylims = [1e-13, 1e-2],
    xticks = [10 * 2^i for i = 0:7],
    yticks = LogTicks(-13:-2),
    xminorgridvisible = false,
    xminorticks = IntervalsBetween(10),
    xscale = log10,
    yscale = log10,
    triangle_bottom = true,
    triangle_top = true,
    triangle_bottom_order = 5,
    triangle_top_order = 5,
    triangle_size = 2.0,
    triangle_shift = 2.0,
    legend_position = (:right, :bottom),
    kwargs...,
)

    # Load data from file 
    set_theme!(Theme(font = font))
    data = Dict(
        dir => CSV.File(
            joinpath(dir, file);
            header = true,
            delim = ' ',
            select = collect(1:7),
            ignorerepeated = true,
        ) for dir in dirs
    )

    # Set up figure parameters
    f = Figure(size = size, fontsize = fontsize, labelfontsize = legendfontsize)
    ax = Axis(
        f[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        xticklabelfont = xticklabelfont,
        yticklabelfont = yticklabelfont,
        xticks = xticks,
        yticks = yticks,
        xscale = xscale,
        yscale = yscale,
        xminorgridvisible = xminorgridvisible,
        xminorticks = xminorticks,
        kwargs...,
    )

    # by default use Makie's automatic axis scaling. Otherwise use custom values.
    if !isnothing(xlims)
        xlims!(ax, xlims)
    end
    if !isnothing(ylims)
        ylims!(ax, ylims)
    end

    # Draw lines for each directory
    for (dir, label, style, color) in zip(dirs, labels, styles, colors)
        if !isnothing(ynorm)
            normalization = data[dir][ynorm]
        else
            normalization = 1.0
        end
        scatterlines!(
            ax,
            data[dir][xkey],
            data[dir][ykey] ./ normalization,
            label = label,
            linestyle = style,
            linewidth = linewidth,
            color = Makie.wong_colors()[color],
        )
    end

    # Make convergence triangle
    if triangle_bottom
        if !isnothing(ynorm)
            normalization = data[dirs[end]][ynorm][end]
        else
            normalization = 1.0
        end
        x0 = data[dirs[end]][xkey][end]
        x1 = x0 * triangle_size
        y0 = (data[dirs[end]][ykey][end] / normalization) / triangle_shift
        y1 = y0 * triangle_size^triangle_bottom_order
        lines!(ax, [x0, x1, x1, x0], [y0, y0, y1, y0]; color = :black)

        # Place text at centroid
        xm = 10^((log10(x0) + 2 * log10(x1)) / 3)
        ym = 10^((2 * log10(y0) + log10(y1)) / 3)
        text!(
            ax,
            string(triangle_bottom_order, ":1");
            position = (xm, ym),
            align = (:center, :center),
            color = :black,
            fontsize = trianglefontsize,
        )
    end

    if triangle_top
        if !isnothing(ynorm)
            normalization = data[dirs[1]][ynorm][end]
        else
            normalization = 1.0
        end
        x0 = data[dirs[1]][xkey][end]
        x1 = x0 * triangle_size
        y0 = (data[dirs[1]][ykey][end] / normalization) * triangle_shift
        y1 = y0 * triangle_size^triangle_top_order
        lines!(ax, [x0, x1, x0, x0], [y0, y1, y1, y0]; color = :black)

        # Place text at centroid (actually shift a bit to be prettier)
        xm = 10^((2.1 * log10(x0) + 0.9 * log10(x1)) / 3)
        ym = 10^((0.9 * log10(y0) + 2.1 * log10(y1)) / 3)
        text!(
            ax,
            string(triangle_top_order, ":1");
            position = (xm, ym),
            align = (:center, :center),
            color = :black,
            fontsize = trianglefontsize,
        )
    end

    if legend
        axislegend(
            ax;
            position = legend_position,
            font = legendfont,
            labelsize = legendfontsize,
        )
    end

    if isnothing(plot_name)
        plot_name = "convergence.pdf"
    end
    save(joinpath(plots_dir, plot_name), f)
end

function plot_evolution(
    dirs = joinpath(
        RESULTS_DIR,
        string(
            Dates.format(today(), dateformat"yyyymmdd"),
            "_steady_barotropic_instability",
        ),
        "N7M4",
    );
    plots_dir = PLOTS_DIR,
    plot_name = nothing,
    labels = [LaTeXString("EC"), LaTeXString("ES")],
    styles = [:dash, :solid],
    colors = 1:length(labels),
    file = "analysis.dat",
    xkey = "time",
    ykey = "entropy",
    ynorm = 1.0,
    xlabel = L"Time $t$ (days)",
    ylabel = LaTeXString("Normalized entropy change"),
    x_in_days = true,
    relative = true,
    font = "CMU Serif",
    size = (350, 350),
    fontsize = 15,
    legendfontsize = 15,
    xticklabelfont = "CMU Serif",
    yticklabelfont = "CMU Serif",
    legendfont = "CMU Serif",
    xlims = nothing,
    ylims = nothing,
    linewidth = 2,
    legend = true,
    legend_position = (:right, :bottom),
    exponent_text = nothing,
    vlinepositions = nothing,
    vlinewidth = 1,
    vlinelabelfont = "CMU Serif",
    vlinelabelsize = 12,
    vlinealigns = [(:right, :bottom), (:right, :bottom)],
    vlineoffsets = [(-2, 0), (-2, 0)],
    vlinelabels = [
        string(LaTeXString("Standard"), "\n", LaTeXString("crash")),
        string(LaTeXString("EC"), "\n", LaTeXString("crash")),
    ],
    vlinecolors = [:black, :black],
    kwargs...,
)

    # Load data from file 
    set_theme!(Theme(font = font))
    data = Dict(
        dir => CSV.File(
            joinpath(dir, file);
            header = true,
            delim = ' ',
            #select = collect(1:18),
            ignorerepeated = true,
        ) for dir in dirs
    )

    # Set up figure parameters
    f = Figure(size = size, fontsize = fontsize)
    ax = Axis(
        f[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        xticklabelfont = xticklabelfont,
        yticklabelfont = yticklabelfont,
        kwargs...,
    )

    # by default use Makie's automatic axis scaling. Otherwise use custom values.
    if !isnothing(xlims)
        xlims!(ax, xlims)
    end
    if !isnothing(ylims)
        ylims!(ax, ylims)
    end

    # Draw lines for each directory
    for (dir, label, style, color) in
        zip(reverse(dirs), reverse(labels), reverse(styles), reverse(colors))
        println("dir = ", dir)
        if x_in_days
            xvalues = data[dir][xkey] / SECONDS_PER_DAY
        else
            xvalues = data[dir][xkey]
        end
        if relative
            yvalues = ((data[dir][ykey] .- data[dir][ykey][1]) / data[dir][ykey][1]) / ynorm
        else
            yvalues = data[dir][ykey] / ynorm
        end
        lines!(
            ax,
            xvalues,
            yvalues,
            label = label,
            linestyle = style,
            linewidth = linewidth,
            color = Makie.Cycled(color),
        )
    end

    if !isnothing(exponent_text)
        Label(f[1, 1, Top()], halign = :left, exponent_text)
    end

    if !isnothing(vlinepositions)
        for i in eachindex(vlinepositions)
            vlines!(ax, [vlinepositions[i]], color = vlinecolors[i], linewidth = vlinewidth)
        end
        if !isnothing(vlinelabels) && !isnothing(ylims)
            for i in eachindex(vlinepositions)
                text!(
                    ax,
                    vlinepositions[i],
                    ylims[1];
                    text = vlinelabels[i],
                    color = vlinecolors[i],
                    align = vlinealigns[i],
                    offset = vlineoffsets[i],
                    font = vlinelabelfont,
                    fontsize = vlinelabelsize,
                )
            end
        end
    end

    if legend
        axislegend(
            ax,
            reverse(ax.scene.plots[1:length(labels)]),
            labels;
            position = legend_position,
            font = legendfont,
            labelsize = legendfontsize,
        )
    end

    if isnothing(plot_name)
        plot_name = string(ykey, "_evolution.pdf")
    end
    save(joinpath(plots_dir, plot_name), f)
end

@inline initial_condition_well_balanced(
    x,
    t,
    equations::TrixiAtmo.AbstractCovariantShallowWaterEquations2D,
) = SVector(5960.0, 0.0, 0.0, 0.0)

@inline function initial_condition_steady_barotropic_instability(x, t, equations)
    RealT = eltype(x)
    a = sqrt(x[1]^2 + x[2]^2 + x[3]^2)  # radius of the sphere
    lon, lat = atan(x[2], x[1]), asin(x[3] / a)

    # compute zonal and meridional velocity components
    u_0 = 80.0f0
    lat_0 = convert(RealT, π / 7)
    lat_1 = convert(RealT, π / 2) - lat_0
    vlon = TrixiAtmo.galewsky_velocity(lat, u_0, lat_0, lat_1)
    vlat = zero(eltype(x))

    # numerically integrate (here we use the QuadGK package) to get height
    galewsky_integral, _ = TrixiAtmo.quadgk(
        latp -> TrixiAtmo.galewsky_integrand(latp, u_0, lat_0, lat_1, a),
        convert(RealT, π / 2),
        lat,
    )
    h = 10158.0f0 - a / EARTH_GRAVITATIONAL_ACCELERATION * galewsky_integral

    # Convert primitive variables from spherical coordinates to the chosen global 
    # coordinate system, which depends on the equation type
    return TrixiAtmo.spherical2global(
        SVector(h, vlon, vlat, zero(RealT), zero(RealT)),
        x,
        equations,
    )
end

# Specialize the L2 and Linf error calculation
function Trixi.calc_error_norms(
    func,
    u,
    t,
    analyzer,
    mesh::P4estMesh{2},
    equations::TrixiAtmo.AbstractCovariantShallowWaterEquations2D,
    initial_condition,
    dg::DGSEM,
    cache,
    cache_analysis,
)
    (; weights) = dg.basis
    (; node_coordinates) = cache.elements
    (; aux_node_vars) = cache.auxiliary_variables

    # Set up data structures
    l2_error = zero(
        func(
            Trixi.get_node_vars(u, equations, dg, 1, 1, 1),
            TrixiAtmo.get_node_aux_vars(aux_node_vars, equations, dg, 1, 1, 1),
            equations,
        ),
    )
    linf_error = copy(l2_error)
    l2_normalization = copy(l2_error)
    linf_normalization = copy(l2_error)

    # Iterate over all elements for error calculations
    for element in eachelement(dg, cache)

        # Calculate errors at each volume quadrature node
        for j in eachnode(dg), i in eachnode(dg)
            x_node = Trixi.get_node_coords(node_coordinates, equations, dg, i, j, element)

            # Convert exact solution into contravariant components using geometric
            # information stored in aux vars
            aux_node =
                TrixiAtmo.get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
            u_exact = initial_condition(x_node, t, aux_node, equations)

            # Compute the difference as usual
            func_exact = func(u_exact, aux_node, equations)
            u_numerical = Trixi.get_node_vars(u, equations, dg, i, j, element)
            diff = func(u_numerical, aux_node, equations) - func_exact

            # For the L2 error, integrate with respect to area element stored in aux vars 
            J = TrixiAtmo.area_element(aux_node, equations)
            l2_error += diff .^ 2 * (weights[i] * weights[j] * J)

            # Compute Linf error as usual
            linf_error = @. max(linf_error, abs(diff))

            # Compute normalization
            linf_normalization = @. max(linf_normalization, abs(func_exact))
            l2_normalization += func_exact .^ 2 * (weights[i] * weights[j] * J)
        end
    end

    # Normalize both errors 
    l2_error = @. sqrt(l2_error / l2_normalization)
    linf_error = @. linf_error / linf_normalization

    return l2_error, linf_error
end

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
        "elixirs/elixir_spherical_shallow_water.jl",
        1,
        polydeg = 7,
        initial_condition = initial_condition_rossby_haurwitz,
        auxiliary_field = nothing,
        surface_flux = surface_flux_ec,
        initial_cells_per_dimension = 8,
        identifier = "_ec_N7",
        interval = 50,
        tspan = (0.0, 28.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 280,
    )

    run_driver(
        "elixirs/elixir_spherical_shallow_water.jl",
        1,
        polydeg = 7,
        initial_condition = initial_condition_rossby_haurwitz,
        auxiliary_field = nothing,
        surface_flux = surface_flux_es,
        initial_cells_per_dimension = 8,
        identifier = "_es_N7",
        interval = 50,
        tspan = (0.0, 28.0 * SECONDS_PER_DAY),
        cfl = 0.1,
        n_saves = 280,
    )
end

function plot_unsteady_solid_body_rotation()
    # Figure 2a
    plot_convergence(
        [
            "../results/20250713_unsteady_solid_body_rotation_ec_N3/",
            "../results/20250713_unsteady_solid_body_rotation_es_N3/",
        ],
        plot_name = "unsteady_solid_body_rotation_convergence_N3.pdf",
        triangle_bottom_order = 4,
        triangle_top_order = 3,
    )
    # Figure 2b
    plot_convergence(
        [
            "../results/20250713_unsteady_solid_body_rotation_ec_N4/",
            "../results/20250713_unsteady_solid_body_rotation_es_N4/",
        ],
        plot_name = "unsteady_solid_body_rotation_convergence_N4.pdf",
        triangle_bottom_order = 5,
        triangle_top_order = 5,
    )

    # Figure 2c
    plot_convergence(
        [
            "../results/20250713_unsteady_solid_body_rotation_ec_p_refine/",
            "../results/20250714_unsteady_solid_body_rotation_es_p_refine/",
        ],
        triangle_bottom = false,
        triangle_top = false,
        plot_name = "unsteady_solid_body_rotation_p_refine_M4.pdf",
        xkey = "N",
        xlabel = L"Polynomial degree $N$",
        xticks = collect(2:10),
        xscale = identity,
        legend_position = (:left, :bottom),
    )
end

function plot_isolated_mountain()
    # Figure 3a
    plot_evolution(
        [
            "../results/20250713_well_balanced_ec/N3M20",
            "../results/20250713_well_balanced_es/N3M20",
        ],
        plot_name = "well_balanced_l2_h_evolution_N3M20.pdf",
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
            "../results/20250713_isolated_mountain_ec/N3M20",
            "../results/20250713_isolated_mountain_es/N3M20",
        ],
        plot_name = "isolated_mountain_mass_evolution_N3M20.pdf",
        legend_position = (:left, :top),
        xlims = [0, 15],
        xticks = [0, 5, 10, 15],
        ylims= [-0.85,4.5],
        ykey = "mass",
        ylabel = LaTeXString("Normalized mass change"),
        ynorm = 1e-14,
        exponent_text = L"\times 10^{-14}",
    )

    # Figure 3c
    plot_evolution(
        [
            "../results/20250713_isolated_mountain_ec/N3M20",
            "../results/20250713_isolated_mountain_es/N3M20",
        ],
        plot_name = "isolated_mountain_entropy_evolution_N3M20.pdf",
        legend_position = (:left, :bottom),
        xlims = [0, 15],
        xticks = [0, 5, 10, 15],
        ylims = [-8.25, 0.25],
        ykey = "entropy",
        ynorm = 1e-8,
        exponent_text = L"\times 10^{-8}",
    )
end

function plot_barotropic_instability()
    # Figure 4a
    plot_evolution(
        [
            "../results/20250525_steady_barotropic_instability_ec_M16/N3M16",
            "../results/20250525_steady_barotropic_instability_es_M16/N3M16",
        ],
        plot_name = "steady_barotropic_instability_l2_h_evolution_N3M16.pdf",
        ykey = "l2_h",
        legend_position = (:left, :top),
        relative = false,
        ynorm = 1e-3,
        exponent_text = L"\times 10^{-3}",
        ylabel = L"Normalized $L^2$ height error",
        xticks = [0, 3, 6, 9, 12],
        xlims = [0, 12],
        ylims = [-1, 12],
    )

    # Figure 4b
    plot_evolution(
        [
            "../results/20250520_steady_barotropic_instability_ec/N3M32",
            "../results/20250520_steady_barotropic_instability_es/N3M32",
        ],
        plot_name = "steady_barotropic_instability_l2_h_evolution_N3M32.pdf",
        ykey = "l2_h",
        legend_position = (:left, :top),
        relative = false,
        ynorm = 1e-3,
        exponent_text = L"\times 10^{-3}",
        ylabel = L"Normalized $L^2$ height error",
        xticks = [0, 3, 6, 9, 12],
        xlims = [0, 12],
        ylims = [-1, 12],
    )

    # Figure 4c
    plot_evolution(
        [
            "../results/20250520_steady_barotropic_instability_ec/N3M64",
            "../results/20250520_steady_barotropic_instability_es/N3M64",
        ],
        plot_name = "steady_barotropic_instability_l2_h_evolution_N3M64.pdf",
        ykey = "l2_h",
        legend_position = (:left, :top),
        relative = false,
        ynorm = 1e-3,
        exponent_text = L"\times 10^{-3}",
        ylabel = L"Normalized $L^2$ height error",
        xticks = [0, 3, 6, 9, 12],
        xlims = [0, 12],
        ylims = [-1, 12],
    )


    # Figure 4d
    plot_evolution(
        [
            "../results/20250525_steady_barotropic_instability_ec_M16/N3M16",
            "../results/20250525_steady_barotropic_instability_es_M16/N3M16",
        ],
        plot_name = "steady_barotropic_instability_entropy_evolution_N3M16.pdf",
        legend_position = (:left, :bottom),
        ykey = "entropy",
        ynorm = 1e-5,
        exponent_text = L"\times 10^{-8}",
        xticks = [0, 3, 6, 9, 12],
        xlims = [0, 12],
        ylims = [-5, 0.25],
    )

    # Figure 4e
    plot_evolution(
        [
            "../results/20250520_steady_barotropic_instability_ec/N3M32",
            "../results/20250520_steady_barotropic_instability_es/N3M32",
        ],
        plot_name = "steady_barotropic_instability_entropy_evolution_N3M32.pdf",
        legend_position = (:left, :bottom),
        ykey = "entropy",
        ynorm = 1e-5,
        exponent_text = L"\times 10^{-8}",
        xticks = [0, 3, 6, 9, 12],
        xlims = [0, 12],
        ylims = [-5, 0.25],
    )

    # Figure 4f
    plot_evolution(
        [
            "../results/20250520_steady_barotropic_instability_ec/N3M64",
            "../results/20250520_steady_barotropic_instability_es/N3M64",
        ],
        plot_name = "steady_barotropic_instability_entropy_evolution_N3M64.pdf",
        legend_position = (:left, :bottom),
        ykey = "entropy",
        ynorm = 1e-5,
        exponent_text = L"\times 10^{-8}",
        xticks = [0, 3, 6, 9, 12],
        xlims = [0, 12],
        ylims = [-5, 0.25],
    )

end

function plot_rossby_haurwitz()
    # Figure 7a
    plot_evolution(
        [
            "../results/20250526_rossby_haurwitz_ec_28days/N3M16/",
            "../results/20250526_rossby_haurwitz_es_28days/N3M16/",
            "../results/20250520_rossby_haurwitz_standard_21days/N3M16/",
        ],
        plot_name = "rossby_haurwitz_entropy_evolution_N3M16.pdf",
        labels = [LaTeXString("EC"), LaTeXString("ES"), LaTeXString("Standard")],
        styles = [:dash, :solid, :dot],
        ykey = "entropy",
        ynorm = 1e-6,
        exponent_text = L"\times 10^{-6}",
        xticks = [0, 7, 14, 21, 28],
        legend_position = (:right, :top),
        xlims = [0, 28],
        ylims = [-5, 5],
        vlinepositions = [18.71370, 26.03088], # crash times
        size = (500, 350),
    )

    # Figure 7b
    plot_evolution(
        [
            "../results/20250520_rossby_haurwitz_ec_21days_N7/N7M8/",
            "../results/20250526_rossby_haurwitz_es_28day_N7/N7M8/",
            "../results/20250520_rossby_haurwitz_standard_21days_N7/N7M8/",
        ],
        plot_name = "rossby_haurwitz_entropy_evolution_N7M8.pdf",
        labels = [LaTeXString("EC"), LaTeXString("ES"), LaTeXString("Standard")],
        styles = [:dash, :solid, :dot],
        ykey = "entropy",
        ynorm = 1e-6,
        exponent_text = L"\times 10^{-6}",
        xticks = [0, 7, 14, 21, 28],
        legend_position = (:right, :top),
        xlims = [0, 28],
        ylims = [-5, 5],
        vlinepositions = [10.91976, 13.92261], # crash times
        size = (500, 350),
    )
end

end # module SphericalShallowWater
