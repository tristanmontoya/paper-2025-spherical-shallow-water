module SphericalShallowWater

using Trixi, TrixiAtmo, Trixi2Vtk
using CairoMakie, LaTeXStrings, Dates, Printf, CSV

export EXAMPLES_DIR, RESULTS_DIR
export run_driver, plot_convergence

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

# plot_convergence(["../results/20250430_unsteady_solid_body_rotation_N4_ec/", "../results/20250430_unsteady_solid_body_rotation_N4_es/"])
function plot_convergence(
    dirs = joinpath(
        RESULTS_DIR,
        string(
            Dates.format(today(), dateformat"yyyymmdd"),
            "_unsteady_solid_body_rotation",
        ),
    );
    labels = ["Entropy conservative", "Entropy stable"],
    file = "analysis.dat",
    xlabel = LaTeXString("Nominal resolution (km)"),
    ylabel = L"Normalized $L^2$ height error",
    font = "CMU Serif",
    size = (400, 300),
    fontsize = 12,
    xticklabelfont = "CMU Serif",
    yticklabelfont = "CMU Serif",
    xlims = nothing,
    ylims = nothing,
    xticks = LogTicks(0:6),
    yticks = LogTicks(-12:-4),
    xscale = log10,
    yscale = log10,
    triangle = true,
    triangle_order = 5,
    triangle_size = 2.0,
    triangle_shift = 0.5,
    kwargs...,
)

    # Load data from file 
    set_theme!(Theme(font = font))
    data = Dict(
        dir => CSV.File(
            joinpath(dir, file);
            header = true,
            delim = ' ',
            select = [3, 4],
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
        kwargs...,
    )

    # by default use Makie's automatic axis scaling. Otherwise use custom values.
    if !isnothing(xlims) || !isnothing(ylims)
        xlims!(ax, xlims)
        ylims!(ax, ylims)
    end

    # Draw lines for each directory
    for (dir, label) in zip(dirs, labels)
        scatterlines!(
            ax,
            data[dir]["resolution_km"],
            data[dir]["l2_height_normalized"],
            label = LaTeXString(label),
        )
    end

    # Make convergence triangle
    if triangle
        x0 = data[dirs[end]]["resolution_km"][end]
        x1 = x0 * triangle_size
        y0 = data[dirs[end]]["l2_height_normalized"][end] * triangle_shift
        y1 = y0 * triangle_size^triangle_order
        lines!(ax, [x0, x1, x1, x0], [y0, y0, y1, y0]; color = :black)

        # Place text at centroid
        xm = 10^((log10(x0) + 2 * log10(x1)) / 3) 
        ym = 10^((2 * log10(y0) + log10(y1)) / 3)
        text!(
            ax,
            string(triangle_order, ":1");
            position = (xm, ym),
            align = (:center, :center),
            color = :black,
        )
    end

    axislegend(ax, position = (:left, :top))
    save(joinpath(first(dirs), "convergence.pdf"), f)
end

function slope_triangle!(
    ax;
    order::Real,
    anchor::Tuple{<:Real,<:Real},
    dx::Real = 2,
    color = :black,
    label::Bool = true,
    labelfontsize::Real = 14,
)

    x0, y0 = anchor                 # lower-left corner in *data* units
    x1 = x0 * dx                # horizontal length  = Δx
    y2 = y0 * dx^(order)       # vertical length    = (Δx)ᵖ   (because log–log)

    println(x0, " ", x1, " ", y0, " ", y2)
    # draw the outline (transparent fill keeps the underlying grid visible)
    poly!(
        ax,
        Point2f[(x0, y0), (x1, y0), (x1, y2)];
        color = :red,
        strokecolor = :black,
        strokewidth = 1.5,
    )
    # :contentReference[oaicite:0]{index=0}
    #=
    # optional “p” label at triangle centre (geometric mean works well)
    if label
        xm = sqrt(x0 * x1)
        ym = sqrt(y0 * y2)
        text!(
            ax,
            string(order);
            position = (xm, ym),
            align = (:center, :center),
            color = color,
            fontsize = labelfontsize,
        )          # :contentReference[oaicite:1]{index=1}
    end
    return nothing
    =#
end

end # module SphericalShallowWater
