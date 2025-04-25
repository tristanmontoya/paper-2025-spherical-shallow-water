module SphericalShallowWater

using Trixi, TrixiAtmo
using Dates

export EXAMPLES_DIR, RESULTS_DIR
export atmo_convergence_test

const EXAMPLES_DIR = TrixiAtmo.examples_dir()
const RESULTS_DIR = joinpath(dirname(dirname(@__DIR__)), "results")

function atmo_convergence_test(
    mod::Module,
    elixir::AbstractString,
    iterations,
    RealT = Float64;
    results_dir = RESULTS_DIR,
    polydeg = 3,
    initial_condition = initial_condition_unsteady_solid_body_rotation,
    identifier = Dates.format(today(), dateformat"yyyymmdd"),
    kwargs...,
)
    @assert(
        iterations > 1,
        "Number of iterations must be bigger than 1 for a convergence analysis"
    )

    # Types of errors to be calculated
    errors = Dict(:l2 => RealT[], :linf => RealT[])

    initial_resolution = Trixi.extract_initial_resolution(elixir, kwargs)

    ic_name = replace(string(initial_condition), r"^initial_condition_" => "")
    # run simulations and extract errors
    for iter = 1:iterations
        println("Running convtest iteration ", iter, "/", iterations)
        cells_per_dimension = initial_resolution .* 2^(iter - 1)
        trixi_include(
            mod,
            elixir;
            kwargs...,
            output_dir = joinpath(
                results_dir,
                string(ic_name, "_", identifier),
                string("N", polydeg, "M", cells_per_dimension[1]),
            ),
            polydeg = polydeg,
            cells_per_dimension = cells_per_dimension,
        )

        l2_error, linf_error = mod.analysis_callback(mod.sol)

        # collect errors as one vector to reshape later
        append!(errors[:l2], l2_error)
        append!(errors[:linf], linf_error)

        println("\n\n")
        println("#"^100)
    end

    # Use raw error values to compute EOC
    Trixi.analyze_convergence(errors, iterations, mod.semi)
end

end # module SphericalShallowWater
