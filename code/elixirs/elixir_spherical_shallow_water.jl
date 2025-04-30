###############################################################################
# Standard elixir for the spherical shallow water equations
###############################################################################

using OrdinaryDiffEq, Trixi, TrixiAtmo

###############################################################################
# Parameters

initial_condition = initial_condition_unsteady_solid_body_rotation
auxiliary_field = bottom_topography_unsteady_solid_body_rotation
polydeg = 3
cells_per_dimension = 2
n_saves = 10
tspan = (0.0, 15.0 * SECONDS_PER_DAY)
output_dir = "out"

###############################################################################
# Custom outputs 
waterheight_squared = (u, aux, equations) -> waterheight(u, equations)^2
Trixi.pretty_form_utf(::typeof(waterheight_squared)) = "∑h²"
Trixi.pretty_form_ascii(::typeof(waterheight_squared)) = "h_squared"

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
    interval = 200,
    save_analysis = true,
    output_directory = output_dir,
    extra_analysis_integrals = (entropy, waterheight_squared),
)

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(
    dt = (tspan[2] - tspan[1]) / n_saves,
    output_directory = output_dir,
    solution_variables = cons2prim_and_vorticity,
)

# The StepsizeCallback handles the re-calculation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl = 0.4)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE 
# solver
callbacks =
    CallbackSet(summary_callback, analysis_callback, save_solution, stepsize_callback)

###############################################################################
# run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed 
# callbacks
sol = solve(
    ode,
    CarpenterKennedy2N54(williamson_condition = false),
    dt = 100.0,
    maxiters=1e7,
    save_everystep = false,
    callback = callbacks,
)

l2_height_error = analysis_callback(sol)[1][1]
l2_height_normalization = sqrt(integrate(waterheight_squared, last(sol.u), semi))
