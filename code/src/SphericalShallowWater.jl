module SphericalShallowWater

using Trixi, TrixiAtmo, Trixi2Vtk
using CairoMakie, LaTeXStrings, Dates, Printf, CSV

export EXAMPLES_DIR, RESULTS_DIR
export run_driver, plot_convergence, plot_evolution
export initial_condition_well_balanced, initial_condition_steady_barotropic_instability

const EXAMPLES_DIR = TrixiAtmo.examples_dir()
const RESULTS_DIR = joinpath(dirname(dirname(@__DIR__)), "results")

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
    fmt1 = Printf.Format("%-4d" * "%-4d" * "%-25.17e"^2 * "missing" * "\n") # no EOC
    fmt2 = Printf.Format("%-4d" * "%-4d" * "%-25.17e"^2 * "%.2f" * "\n")
    headers = ["N   ", "M   ", "resolution_km", "l2_height_normalized", "order"]
    open(joinpath(project_dir, "analysis.dat"), "w") do io
        println(io, string(headers[1:2]..., rpad.(headers[3:end-1], 25)..., headers[end]))
    end

    resolutions = RealT[]
    errors = RealT[]

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

            # Note that by default, both the error and normalization factor have been 
            # scaled by the domain size.
            l2_height_normalized = mod.l2_height_error / mod.l2_height_normalization
            resolution_km = π * EARTH_RADIUS / (M * N * 1000) # scale by 1000 to get km

            append!(resolutions, resolution_km)
            append!(errors, l2_height_normalized)

            open(joinpath(project_dir, "analysis.dat"), "a") do io
                if M == initial_cells_per_dimension # don't compute order for first grid
                    Printf.format(io, fmt1, N, M, resolution_km, l2_height_normalized)
                else
                    Printf.format(
                        io,
                        fmt2,
                        N,
                        M,
                        resolution_km,
                        l2_height_normalized,
                        log(errors[end] / errors[end-1]) / log(1 / 2),
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
    labels = [LaTeXString("Entropy conservative"), LaTeXString("Entropy stable")],
    styles = [:dash, :solid],
    colors = [:black, :black],
    file = "analysis.dat",
    xkey = "resolution_km",
    ykey = "l2_height_normalized",
    xlabel = LaTeXString("Nominal resolution (km)"),
    ylabel = L"Normalized $L^2$ height error",
    font = "CMU Serif",
    size = (350, 350),
    fontsize = 12,
    xticklabelfont = "CMU Serif",
    yticklabelfont = "CMU Serif",
    legendfont = "CMU Serif",
    xlims = nothing,
    ylims = [1e-12, 1e-2],
    xticks = LogTicks(1:4),
    yticks = LogTicks(-12:-2),
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
            select = [1, 2, 3, 4],
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
        xticks = xticks,
        yticks = yticks,
        xscale = xscale,
        yscale = yscale,
        xminorgridvisible = true,
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
        scatterlines!(
            ax,
            data[dir][xkey],
            data[dir][ykey],
            label = label,
            linestyle = style,
            color = color,
        )
    end

    # Make convergence triangle
    if triangle_bottom
        x0 = data[dirs[end]][xkey][end]
        x1 = x0 * triangle_size
        y0 = data[dirs[end]][ykey][end] / triangle_shift
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
        )
    end

    if triangle_top
        x0 = data[dirs[1]][xkey][end]
        x1 = x0 * triangle_size
        y0 = data[dirs[1]][ykey][end] * triangle_shift
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
        )
    end

    axislegend(ax, position = legend_position, font = legendfont)
    save(joinpath(first(dirs), "convergence.pdf"), f)
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
    labels = [LaTeXString("Entropy conservative"), LaTeXString("Entropy stable")],
    styles = [:dash, :solid],
    colors = [:black, :black],
    file = "analysis.dat",
    xkey = "time",
    ykey = "entropy",
    xlabel = LaTeXString("Time (days)"),
    ylabel = LaTeXString("Normalized entropy change"),
    x_in_days = true,
    relative = true,
    font = "CMU Serif",
    size = (350, 350),
    fontsize = 12,
    xticklabelfont = "CMU Serif",
    yticklabelfont = "CMU Serif",
    legendfont = "CMU Serif",
    xlims = nothing,
    ylims = nothing,
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
            select = collect(1:12),
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
    for (dir, label, style, color) in zip(dirs, labels, styles, colors)
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
        lines!(ax, xvalues, yvalues, label = label, linestyle = style, color = color)
    end

    axislegend(ax, position = legend_position, font = legendfont)
    save(joinpath(first(dirs), string(ykey, "_evolution.pdf")), f)
end

@inline initial_condition_well_balanced(x, t, equations) = SVector(5960.0, 0.0, 0.0, 0.0)

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
