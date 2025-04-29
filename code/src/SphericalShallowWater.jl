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
    initial_resolution = 4, # initial value of `cells_per_dimension`
    initial_condition = initial_condition_unsteady_solid_body_rotation,
    date = Dates.format(today(), dateformat"yyyymmdd"),
    identifier = "", # suffix to identify a specific run
    kwargs...,
)

    mod = @__MODULE__

    @assert(
        iterations > 1,
        "Number of iterations must be bigger than 1 for a convergence analysis"
    )

    # Get name for project
    ic_name = replace(string(initial_condition), r"^initial_condition_" => "")
    project_dir = mkpath(joinpath(results_dir, string(date, "_", ic_name, identifier)))
    println("Project directory: ", project_dir)

    # Format top-level output file and write headers
    fmt1 = Printf.Format("%-4d" * "%-4d" * "%-25.17e"^2 * "missing" * "\n") # no EOC
    fmt2 = Printf.Format("%-4d" * "%-4d" * "%-25.17e"^2 * "%.2f" * "\n")
    headers = ["N   ", "M   ", "resolution", "l2_height_normalized", "order"]
    open(joinpath(project_dir, "analysis.dat"), "w") do io
        println(io, string(headers[1:2]..., rpad.(headers[3:end-1], 25)..., headers[end]))
    end

    resolutions = RealT[]
    errors = RealT[]

    cells_per_dimension = initial_resolution .* 2 .^ ((1:iterations) .- 1)

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
            resolution = π * EARTH_RADIUS / (M * N)

            append!(resolutions, resolution)
            append!(errors, l2_height_normalized)

            open(joinpath(project_dir, "analysis.dat"), "a") do io
                if M == initial_resolution # don't compute order for first grid
                    Printf.format(io, fmt1, N, M, resolution, l2_height_normalized)
                else
                    Printf.format(
                        io,
                        fmt2,
                        N,
                        M,
                        resolution,
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
    dir = joinpath(
        RESULTS_DIR,
        string(
            Dates.format(today(), dateformat"yyyymmdd"),
            "_unsteady_solid_body_rotation",
        ),
    );
    file = "analysis.dat",
    xlabel = LaTeXString("Nominal resolution (km)"),
    ylabel = L"Normalized $L^2$ height error",
    font = "CMU Serif",
    xticklabelfont = "CMU Serif",
    yticklabelfont = "CMU Serif",
    xlims = nothing,
    ylims = nothing,
    xticks = LogTicks(3:7),
    yticks = LogTicks(-5:-2),
    xscale = log10,
    yscale = log10,
    kwargs...,
)

    set_theme!(Theme(font = font))
    data = CSV.File(
        joinpath(dir, file);
        header = true,
        delim = ' ',
        select = [3, 4],
        ignorerepeated = true,
    )

    f = Figure(size = (400, 300))
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

    scatterlines!(ax, data["resolution"], data["l2_height_normalized"])
    save(joinpath(dir, "convergence.pdf"), f)
end

end # module SphericalShallowWater
