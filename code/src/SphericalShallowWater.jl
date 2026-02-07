module SphericalShallowWater

using Trixi, TrixiAtmo, Trixi2Vtk
using CairoMakie, LaTeXStrings, Dates, Printf, CSV
using LinearAlgebra: norm

export EXAMPLES_DIR, RESULTS_DIR
export surface_flux_ec, surface_flux_es
export pot_enst
export run_driver, plot_convergence, plot_evolution, calc_norms
export run_unsteady_solid_body_rotation,
    run_isolated_mountain, run_barotropic_instability, run_rossby_haurwitz
export run_timestep_study, plot_timestep_study
export plot_unsteady_solid_body_rotation,
    plot_isolated_mountain, plot_barotropic_instability, plot_rossby_haurwitz
export initial_condition_well_balanced, initial_condition_steady_barotropic_instability

const EXAMPLES_DIR = TrixiAtmo.examples_dir()
const RESULTS_DIR = joinpath(dirname(dirname(@__DIR__)), "results")
const PLOTS_DIR = joinpath(dirname(dirname(@__DIR__)), "plots")

# extra analysis for normalized errors
include("analysis.jl")

# scripts to run drivers
include("run.jl")

# scripts to generate figures
include("plot.jl")

# EC and ES surface fluxes
const surface_flux_ec = (flux_ec, flux_nonconservative_surface_simplified)
const surface_flux_es = (
    FluxPlusDissipation(flux_ec, DissipationLocalLaxFriedrichs()),
    flux_nonconservative_surface_simplified,
)

function run_driver(
    elixir::AbstractString,
    iterations = 1,
    RealT = Float64;
    results_dir = RESULTS_DIR,
    polydeg = 3,
    initial_cells_per_dimension = 4,
    initial_condition = initial_condition_unsteady_solid_body_rotation,
    date = Dates.format(today(), dateformat"yyyymmdd"),
    identifier = "",
    kwargs...,
)

    mod = @__MODULE__

    # Get name for project
    ic_name = replace(string(initial_condition), r"^initial_condition_" => "")
    project_dir = mkpath(joinpath(results_dir, string(date, "_", ic_name, identifier)))
    println("Project directory: ", project_dir)

    # Format top-level output file and write headers
    fmt1 = Printf.Format(
        "%-4d" * "%-4d" * "%-9.5f" * "%-25.17e"^5 * "missing  " * "missing  " * "\n",
    ) # no EOC
    fmt2 = Printf.Format(
        "%-4d" * "%-4d" * "%-9.5f" * "%-25.17e"^5 * "%-9.2f" * "%-9.2f" * "\n",
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
                rpad.(headers[4:(end-2)], 25)...,
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
    polydeg_iter = polydeg isa Integer ? (polydeg:polydeg) : polydeg
    for N in polydeg_iter
        empty!(resolutions)
        empty!(l2_errors)
        empty!(linf_errors)
        empty!(end_times)
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

            # invokelatest used to avoid world age issues
            t_final = Base.invokelatest(getproperty, mod, :t_final)
            l2_depth_error = Base.invokelatest(getproperty, mod, :l2_depth_error)
            linf_depth_error = Base.invokelatest(getproperty, mod, :linf_depth_error)
            l2_height_error = Base.invokelatest(getproperty, mod, :l2_height_error)
            linf_height_error = Base.invokelatest(getproperty, mod, :linf_height_error)

            end_time = t_final / SECONDS_PER_DAY

            push!(end_times, end_time)
            push!(resolutions, resolution_km)
            push!(l2_errors, l2_height_error)
            push!(linf_errors, linf_height_error)

            open(joinpath(project_dir, "analysis.dat"), "a") do io
                if M == first(cells_per_dimension) # don't compute order for first grid
                    Printf.format(
                        io,
                        fmt1,
                        N,
                        M,
                        end_time,
                        resolution_km,
                        l2_depth_error,
                        linf_depth_error,
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
                        l2_depth_error,
                        linf_depth_error,
                        l2_height_error,
                        linf_height_error,
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

function run_timestep_study(
    elixir::AbstractString,
    iterations = 6,
    RealT = Float64;
    results_dir = RESULTS_DIR,
    polydeg = 3,
    cells_per_dimension = 16,
    n_saves = 120,
    initial_condition = initial_condition_steady_barotropic_instability,
    auxiliary_field = nothing,
    surface_flux = surface_flux_ec,
    tspan = (0.0, 12.0 * SECONDS_PER_DAY),
    initial_dt = 320.0,
    date = Dates.format(today(), dateformat"yyyymmdd"),
    identifier = "_timestep_study",
    kwargs...,
)

    mod = @__MODULE__

    # Get name for project
    ic_name = replace(string(initial_condition), r"^initial_condition_" => "")
    project_dir = mkpath(joinpath(results_dir, string(date, "_", ic_name, identifier)))
    println("Project directory: ", project_dir)

    # Format output file and write headers
    fmt1 = Printf.Format(
        "%-4d" *
        "%-4d" *
        "%-25.17e" * # dt
        "%-25.17e" * # mass_initial
        "%-25.17e" * # mass_final
        "%-25.17e" * # mass_error
        "%-25s" *    # mass_order
        "%-25.17e" * # entropy_initial
        "%-25.17e" * # entropy_final
        "%-25.17e" * # entropy_error
        "%-25s" *    # entropy_order
        "\n",
    ) # no EOC for first iteration
    fmt2 = Printf.Format(
        "%-4d" *
        "%-4d" *
        "%-25.17e" * # dt
        "%-25.17e" * # mass_initial
        "%-25.17e" * # mass_final
        "%-25.17e" * # mass_error
        "%-25.2f" *  # mass_order
        "%-25.17e" * # entropy_initial
        "%-25.17e" * # entropy_final
        "%-25.17e" * # entropy_error
        "%-25.2f" *  # entropy_order
        "\n",
    )
    headers = [
        "N   ",
        "M   ",
        "dt                       ",
        "mass_initial             ",
        "mass_final               ",
        "mass_error               ",
        "mass_order               ",
        "entropy_initial          ",
        "entropy_final            ",
        "entropy_error            ",
        "entropy_order            ",
    ]
    open(joinpath(project_dir, "timestep_analysis.dat"), "w") do io
        println(io, string(headers...))
    end

    dt_values = RealT[]
    mass_errors = RealT[]
    entropy_errors = RealT[]

    # Run simulations with decreasing time steps
    for (iter, dt_factor) in enumerate(2.0 .^ ((0:(iterations-1))))
        dt = initial_dt / dt_factor

        # Run simulation with fixed time step (no CFL-based adaptation)
        trixi_include(
            mod,
            elixir;
            kwargs...,
            output_dir = joinpath(
                project_dir,
                string("N", polydeg, "M", cells_per_dimension, "_dt", dt),
            ),
            initial_condition = initial_condition,
            auxiliary_field = auxiliary_field,
            surface_flux = surface_flux,
            polydeg = polydeg,
            cells_per_dimension = cells_per_dimension,
            tspan = tspan,
            dt_initial = dt,
            adapt_timestep = false,
            interval = 50,
            n_saves = n_saves,
        )

        analysis_file = joinpath(
            project_dir,
            string("N", polydeg, "M", cells_per_dimension, "_dt", dt),
            "analysis.dat",
        )

        data = CSV.File(analysis_file; header = true, delim = ' ', ignorerepeated = true)

        # Extract mass and entropy errors
        mass_initial = data.mass[1]
        mass_final = data.mass[end]
        mass_error = abs(mass_final - mass_initial) / abs(mass_initial)
        entropy_initial = data.entropy[1]
        entropy_final = data.entropy[end]
        entropy_error = abs(entropy_final - entropy_initial) / abs(entropy_initial)

        push!(dt_values, dt)
        push!(entropy_errors, entropy_error)
        push!(mass_errors, mass_error)

        # Write results
        open(joinpath(project_dir, "timestep_analysis.dat"), "a") do io
            if iter == 1 # first iteration, no order calculation
                Printf.format(
                    io,
                    fmt1,
                    polydeg,
                    cells_per_dimension,
                    dt,
                    mass_initial,
                    mass_final,
                    mass_error,
                    "missing",
                    entropy_initial,
                    entropy_final,
                    entropy_error,
                    "missing",
                )
            else
                mass_order =
                    log(mass_errors[end] / mass_errors[end-1]) /
                    log(dt_values[end] / dt_values[end-1])
                entropy_order =
                    log(entropy_errors[end] / entropy_errors[end-1]) /
                    log(dt_values[end] / dt_values[end-1])
                Printf.format(
                    io,
                    fmt2,
                    polydeg,
                    cells_per_dimension,
                    dt,
                    mass_initial,
                    mass_final,
                    mass_error,
                    mass_order,
                    entropy_initial,
                    entropy_final,
                    entropy_error,
                    entropy_order,
                )
            end
        end
        println("\n\n")
        println("#"^100)
    end
end

function plot_convergence(
    dirs = [
        joinpath(
            RESULTS_DIR,
            string(
                Dates.format(today(), dateformat"yyyymmdd"),
                "_unsteady_solid_body_rotation",
            ),
        ),
    ];
    plots_dir = PLOTS_DIR,
    plot_name = nothing,
    labels = [LaTeXString("EC"), LaTeXString("ES")],
    styles = [:dash, :solid],
    colors = 1:length(labels),
    file = "analysis.dat",
    select_cols = collect(1:7),
    xkey = "resolution_km",
    ykey = "l2_height_error",
    ykeys = nothing,
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
    triangle_bottom_attach_index = nothing,
    legend_position = (:right, :bottom),
    kwargs...,
)

    mkpath(plots_dir)

    # Load data from file 
    set_theme!(Theme(font = font))
    data = Dict(
        dir => begin
            path = joinpath(dir, file)
            if isnothing(select_cols)
                CSV.File(path; header = true, delim = ' ', ignorerepeated = true)
            else
                CSV.File(
                    path;
                    header = true,
                    delim = ' ',
                    select = select_cols,
                    ignorerepeated = true,
                )
            end
        end for dir in dirs
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

    # Draw lines
    if isnothing(ykeys)
        # Standard mode: one curve per directory
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
    else
        # Multi-y mode: multiple curves from a single directory
        if length(dirs) != 1
            error("plot_convergence: ykeys requires exactly one directory in dirs")
        end
        ykey = ykeys[1] # for triangles
        dir = first(dirs)
        for (yk, label, style, color) in zip(ykeys, labels, styles, colors)
            if !isnothing(ynorm)
                normalization = data[dir][ynorm]
            else
                normalization = 1.0
            end
            scatterlines!(
                ax,
                data[dir][xkey],
                data[dir][yk] ./ normalization,
                label = label,
                linestyle = style,
                linewidth = linewidth,
                color = Makie.wong_colors()[color],
            )
        end
    end

    # Make convergence triangle
    if triangle_bottom
        attach_index = if isnothing(triangle_bottom_attach_index)
            isnothing(ykeys) ? length(dirs) : 1
        else
            triangle_bottom_attach_index
        end

        dir_attach = if isnothing(ykeys)
            if !(1 <= attach_index <= length(dirs))
                error(
                    "plot_convergence: triangle_bottom_attach_index=$(attach_index) out of range. " *
                    "Expected 1:$(length(dirs)).",
                )
            end
            dirs[attach_index]
        else
            if !(1 <= attach_index <= length(ykeys))
                error(
                    "plot_convergence: triangle_bottom_attach_index=$(attach_index) out of range for ykeys. " *
                    "Expected 1:$(length(ykeys)).",
                )
            end
            first(dirs)
        end

        ykey_attach = isnothing(ykeys) ? ykey : ykeys[attach_index]

        if !isnothing(ynorm)
            normalization = data[dir_attach][ynorm][end]
        else
            normalization = 1.0
        end
        x0 = data[dir_attach][xkey][end]
        x1 = x0 * triangle_size
        y0 = (data[dir_attach][ykey_attach][end] / normalization) / triangle_shift
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
    dirs = [
        joinpath(
            RESULTS_DIR,
            string(
                Dates.format(today(), dateformat"yyyymmdd"),
                "_steady_barotropic_instability",
            ),
            "N7M4",
        ),
    ];
    line_order = [2, 1],
    plots_dir = PLOTS_DIR,
    plot_name = nothing,
    plot_absolute = false,
    labels = [LaTeXString("EC"), LaTeXString("ES")],
    styles = [:dot, :solid],
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
    verbose = false,
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
    inset_ylims = nothing,
    inset_xlims = nothing,
    inset_xticks = nothing,
    inset_yticks = nothing,
    inset_halign = 0.025,
    inset_valign = 0.96,
    inset_yaxisposition = :right,
    show_in_inset = [2],
    reverse_foreground_order = false,
    kwargs...,
)

    mkpath(plots_dir)

    # Guard defaults when only a single directory is provided.
    if maximum(line_order) > length(dirs)
        line_order = collect(1:length(dirs))
    end
    if maximum(show_in_inset) > length(dirs)
        show_in_inset = [1]
    end

    # Load data from file 
    set_theme!(Theme(font = font))
    data = Dict(
        dir => CSV.File(
            joinpath(dir, file);
            header = true,
            delim = ' ',
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

    # draw vertical lines for crash times 
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
                    label = nothing,
                )
            end
        end
    end

    # Draw lines for each directory
    ordered_dirs = dirs[line_order]
    ordered_labels = labels[line_order]
    ordered_styles = styles[line_order]
    ordered_colors = colors[line_order]
    nlines = length(ordered_dirs)
    for (i, (dir, label, style, color)) in
        enumerate(zip(ordered_dirs, ordered_labels, ordered_styles, ordered_colors))
        if x_in_days
            xvalues = data[dir][xkey] / SECONDS_PER_DAY
        else
            xvalues = data[dir][xkey]
        end
        yvalues = if relative
            (data[dir][ykey] .- data[dir][ykey][1]) / data[dir][ykey][1]
        else
            data[dir][ykey]
        end

        if plot_absolute
            yvalues = abs.(yvalues)
        end

        yvalues = yvalues / ynorm
        plt = lines!(
            ax,
            xvalues,
            yvalues,
            label = label,
            linestyle = style,
            linewidth = linewidth,
            color = (color isa Integer ? Makie.Cycled(color) : color),
        )
        if reverse_foreground_order
            # Keep legend/order semantics intact (based on insertion order), but flip
            # foreground/background stacking via a z-translation.
            translate!(plt, 0, 0, nlines - i)
        end
        if verbose
            finite_y = filter(isfinite, yvalues)
            if isempty(finite_y)
                println("dir: ", dir, "\nmin/max: (no finite y-values)")
            else
                println(
                    "dir: ",
                    dir,
                    "\nmin: ",
                    minimum(finite_y),
                    ", max: ",
                    maximum(finite_y),
                )
            end
        end
    end
    if verbose
        println("\n")
    end

    if !isnothing(exponent_text)
        Label(f[1, 1, Top()], halign = :left, exponent_text)
    end

    if !isnothing(inset_ylims)
        ax_inset = Axis(
            f[1, 1],
            width = Relative(0.3),
            height = Relative(0.4),
            halign = inset_halign,
            valign = inset_valign,
            xticklabelfont = xticklabelfont,
            yticklabelfont = yticklabelfont,
            yaxisposition = inset_yaxisposition,
            tellheight = false,
            tellwidth = false,
            alignmode = Mixed(
                left = Makie.Protrusion(0),
                right = Makie.Protrusion(0),
                bottom = Makie.Protrusion(0),
                top = Makie.Protrusion(0),
            ),
            xticks = inset_xticks,
            yticks = inset_yticks,
        )
        ylims!(ax_inset, inset_ylims)
        xlims!(ax_inset, inset_xlims)
        inset_dirs = dirs[show_in_inset][line_order[show_in_inset]]
        inset_labels = labels[show_in_inset][line_order[show_in_inset]]
        inset_styles = styles[show_in_inset][line_order[show_in_inset]]
        inset_colors = colors[show_in_inset][line_order[show_in_inset]]
        ninset = length(inset_dirs)
        for (i, (dir, label, style, color)) in
            enumerate(zip(inset_dirs, inset_labels, inset_styles, inset_colors))
            if x_in_days
                xvalues = data[dir][xkey] / SECONDS_PER_DAY
            else
                xvalues = data[dir][xkey]
            end

            if relative
                yvalues = (data[dir][ykey] .- data[dir][ykey][1]) / data[dir][ykey][1]
            else
                yvalues = data[dir][ykey]
            end

            if plot_absolute
                yvalues = abs.(yvalues)
            end

            plt_inset = lines!(
                ax_inset,
                xvalues,
                yvalues,
                label = label,
                linestyle = style,
                linewidth = linewidth,
                color = (color isa Integer ? Makie.Cycled(color) : color),
            )
            if reverse_foreground_order
                translate!(plt_inset, 0, 0, ninset - i)
            end
            translate!(ax_inset.blockscene, 0, 0, 150)
        end
    end

    if legend
        legend_plots = if isnothing(vlinepositions)
            ax.scene.plots
        else
            ax.scene.plots[(length(vlinepositions)*2+1):end]
        end
        leg = axislegend(
            ax,
            legend_plots[line_order],
            labels;
            position = legend_position,
            font = legendfont,
            labelsize = legendfontsize,
        )
        # Keep legend above all curves regardless of any z-order tweaks to the lines.
        translate!(leg.blockscene, 0, 0, 1000)
    end
    resize_to_layout!(f)
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

end # module SphericalShallowWater
