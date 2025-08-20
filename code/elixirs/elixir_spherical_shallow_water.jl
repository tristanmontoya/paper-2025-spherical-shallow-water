###############################################################################
# Standard elixir for the covariant shallow water equations on the sphere using
# entropy-conservative entropy-stable schemes
###############################################################################

using OrdinaryDiffEq, Trixi, TrixiAtmo

###############################################################################
# Parameters

initial_condition = initial_condition_unsteady_solid_body_rotation
auxiliary_field = bottom_topography_unsteady_solid_body_rotation
polydeg = 3
cells_per_dimension = 2
n_saves = 50
tspan = (0.0, 5.0 * SECONDS_PER_DAY)
output_dir = "out"

###############################################################################
# Custom outputs
mass = (u, aux, equations) -> waterheight(u, equations)
Trixi.pretty_form_utf(::typeof(mass)) = "mass"
Trixi.pretty_form_ascii(::typeof(mass)) = "mass"

###############################################################################
# Spatial discretization

mesh = P4estMeshCubedSphere2D(
    cells_per_dimension,
    EARTH_RADIUS,
    polydeg = polydeg,
    element_local_mapping = true,
)

equations = SplitCovariantShallowWaterEquations2D(
    EARTH_GRAVITATIONAL_ACCELERATION,
    EARTH_ROTATION_RATE,
    global_coordinate_system = GlobalCartesianCoordinates(),
)

# Use entropy-conservative two-point fluxes for volume and surface terms, with the surface 
# flux simplified due to the continuous bottom topography
volume_flux = (flux_ec, flux_nonconservative_ec)
surface_flux = (flux_ec, flux_nonconservative_surface_simplified)

# Create DG solver with polynomial degree = polydeg
solver = DGSEM(
    polydeg = polydeg,
    surface_flux = surface_flux,
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux),
)

# Transform the initial condition to the proper set of conservative variables
initial_condition_transformed = transform_initial_condition(initial_condition, equations)

# A semidiscretization collects data structures and functions for the spatial 
# discretization. Here, we pass in the additional keyword argument "auxiliary_field" to 
# specify the bottom topography.
semi = SemidiscretizationHyperbolic(
    mesh,
    equations,
    initial_condition_transformed,
    solver,
    source_terms = source_terms_geometric_coriolis,
    auxiliary_field = auxiliary_field,
)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0 to T
ode = semidiscretize(semi, tspan)

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation 
# setup and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the
# results. Note that entropy should be conserved at the semi-discrete level.
analysis_callback = AnalysisCallback(
    semi,
    interval = 50,
    save_analysis = true,
    output_directory = output_dir,
    extra_analysis_errors = (:l2_error_primitive, :linf_error_primitive),
    extra_analysis_integrals = (mass, entropy, pot_enst),
)

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(
    dt = (tspan[2] - tspan[1]) / n_saves,
    output_directory = output_dir,
    solution_variables = cons2prim_and_vorticity,
)

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl = 0.1)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE 
# solver
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
    integrator.u,
    integrator.t,
    analysis_callback.affect!.analyzer,
    semi,
    analysis_callback.affect!.cache,
)
t_final = integrator.t
l2_depth_error, linf_depth_error = l2_error[1], linf_error[1]
l2_height_error, linf_height_error = l2_error_prim[1], linf_error_prim[1]
