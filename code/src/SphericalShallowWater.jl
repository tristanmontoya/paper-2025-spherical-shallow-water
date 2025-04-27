module SphericalShallowWater

using Trixi, TrixiAtmo
using Dates, Printf

export EXAMPLES_DIR, RESULTS_DIR
export driver_convergence_test

const EXAMPLES_DIR = TrixiAtmo.examples_dir()
const RESULTS_DIR = joinpath(dirname(dirname(@__DIR__)), "results")

function driver_convergence_test(
    elixir::AbstractString,
    iterations = 1,
    RealT = Float64;
    results_dir = RESULTS_DIR,
    polydeg = 3:4,
    initial_resolution = 4,
    initial_condition = initial_condition_unsteady_solid_body_rotation,
    date = Dates.format(today(), dateformat"yyyymmdd"),
    identifier = "",
    kwargs...,
)

    mod = @__MODULE__

    @assert(
        iterations > 1,
        "Number of iterations must be bigger than 1 for a convergence analysis"
    )

    # Get name for project
    ic_name = replace(string(initial_condition), r"^initial_condition_" => "")
    project_dir =
        mkpath(joinpath(results_dir, string(date, "_", ic_name, "_conv", identifier)))

    # Format output file and write headers
    fmt = Printf.Format("%-4d" * "%-4d" * "%-25.17e"^3 * "\n")
    headers = ["N   ", "M   ", "Resolution (km)", "l2 height error", "linf height error"]
    open(joinpath(project_dir, "analysis.dat"), "w") do io
        println(io, string(headers[1:2]..., rpad.(headers[3:end], 25)..., " "))
    end

    # Types of errors to be calculated
    resolutions = RealT[]
    errors = Dict(:l2 => RealT[], :linf => RealT[])

    cells_per_dimension = initial_resolution .* 2 .^ ((1:iterations) .-1)

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

            l2_error, linf_error = mod.analysis_callback(mod.sol)
            resolution = π * EARTH_RADIUS / (M * N)

            # collect errors as one vector to reshape later
            append!(resolutions, resolution)
            append!(errors[:l2], l2_error)
            append!(errors[:linf], linf_error)

            open(joinpath(project_dir, "analysis.dat"), "a") do io
                Printf.format(
                    io,
                    fmt,
                    N,
                    M,
                    resolution,
                    l2_error[1],
                    linf_error[1],
                )
            end

            println("\n\n")
            println("#"^100)
        end
    end

end

end # module SphericalShallowWater
