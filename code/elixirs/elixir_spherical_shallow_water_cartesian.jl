
using OrdinaryDiffEq
using Trixi
using TrixiAtmo

###############################################################################
# Entropy conservation for the spherical shallow water equations in Cartesian
# form obtained through the projection of the momentum onto the divergence-free
# tangential contravariant vectors

###############################################################################
# Parameters

initial_condition = initial_condition_unsteady_solid_body_rotation
polydeg = 3
cells_per_dimension = 2
n_saves = 50
tspan = (0.0, 5.0 * SECONDS_PER_DAY)
output_dir = "out"
###############################################################################
# Spatial discretization
equations = ShallowWaterEquations3D(
    gravity = EARTH_GRAVITATIONAL_ACCELERATION,
    rotation_rate = EARTH_ROTATION_RATE,
)

# Create DG solver with polynomial degree = 3 and Wintemeyer et al.'s flux as surface flux
volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)

# For provably entropy-stable surface fluxes, use
# surface_flux = (FluxPlusDissipation(flux_wintermeyer_etal, DissipationLaxFriedrichsEntropyVariables()), 
#                 flux_nonconservative_wintermeyer_etal)

solver = DGSEM(
    polydeg = polydeg,
    surface_flux = surface_flux,
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux),
)

# Transform the initial condition to the proper set of conservative variables
initial_condition_transformed = transform_initial_condition(initial_condition, equations)

mesh = P4estMeshCubedSphere2D(
    cells_per_dimension[1],
    EARTH_RADIUS,
    polydeg = polydeg,
    initial_refinement_level = 0,
)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(
    mesh,
    equations,
    initial_condition_transformed,
    solver,
    metric_terms = MetricTermsInvariantCurl(),
    source_terms = source_terms_coriolis_lagrange_multiplier,
)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0.0 to π
ode = semidiscretize(semi, tspan)

# Clean the initial condition
for element in eachelement(solver, semi.cache)
    for j in eachnode(solver), i in eachnode(solver)
        u0 = Trixi.wrap_array(ode.u0, semi)

        contravariant_normal_vector = Trixi.get_contravariant_vector(
            3,
            semi.cache.elements.contravariant_vectors,
            i,
            j,
            element,
        )
        clean_solution_lagrange_multiplier!(
            u0[:, i, j, element],
            equations,
            contravariant_normal_vector,
        )
    end
end

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(
    semi,
    interval = 10,
    save_analysis = true,
    extra_analysis_errors = (:l2_error_primitive, :linf_error_primitive),
    extra_analysis_integrals = (waterheight, energy_total),
    analysis_polydeg = polydeg,
)

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(
    dt = (tspan[2] - tspan[1]) / n_saves,
    output_directory = output_dir,
    solution_variables = cons2prim_and_vorticity,
)

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl = 0.7)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks =
    CallbackSet(summary_callback, analysis_callback, save_solution, stepsize_callback)

###############################################################################
# run the simulation

# Set up integrator
integrator = init(
    ode,
    CarpenterKennedy2N54(williamson_condition = false, thread = Trixi.True()),
    dt = 100.0,
    maxiters = 1e8,
    save_everystep = false,
    callback = callbacks,
)

# get initial residual and save to file for plotting
du = similar(ode.u0)
Trixi.rhs!(du, ode.u0, semi, tspan[1])
save_residual =
    SaveSolutionCallback(output_directory = output_dir, solution_variables = cons2cons)
Trixi.save_solution_file(semi, du, save_residual.affect!, integrator, system = "residual")

# Try to solve; if it fails, detect crash gracefully. This may not work with threading.
try
    solve!(integrator)
catch error
    println("Crashed at t = ", integrator.t)
end

# Calculate error norms at end of simulation
l2_error, linf_error = analysis_callback(integrator.sol)
l2_error_prim, linf_error_prim = Trixi.calc_error_norms(
    cons2prim,
    integrator.t,
    analysis_callback.analyzer,
    semi,
    cache_analysis,
)
t_final = integrator.t
l2_depth_error, linf_depth_error = l2_error[1], linf_error[1]
l2_height_error, linf_height_error = l2_error_prim[1], linf_error_prim[1]
