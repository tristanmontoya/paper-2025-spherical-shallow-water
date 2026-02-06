function pot_enst end
Trixi.pretty_form_utf(::typeof(pot_enst)) = "pot_enst"
Trixi.pretty_form_ascii(::typeof(pot_enst)) = "pot_enst"

# Specialize the L2 and Linf error calculation, since Trixi.jl does not normalize the same 
# way as Williamson et al. (1992) and other geophysical fluid dynamics papers
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

# Potential enstrophy
function Trixi.analyze(
    ::typeof(pot_enst),
    du,
    u,
    t,
    mesh::P4estMesh{2},
    equations::TrixiAtmo.AbstractCovariantShallowWaterEquations2D,
    dg::DGSEM,
    cache,
)
    (; aux_node_vars,) = cache.auxiliary_variables
    (; node_coordinates,) = cache.elements

    Trixi.integrate_via_indices(
        u,
        mesh,
        equations,
        dg,
        cache,
        du,
    ) do u, i, j, element, equations, dg, du_node
        x_node = Trixi.get_node_coords(node_coordinates, equations, dg, i, j, element)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)

        h = Trixi.waterheight(u_node, equations)
        f = 2 * equations.rotation_rate * x_node[3] / norm(x_node)  # 2Ωsinθ
        zeta = TrixiAtmo.calc_vorticity_node(u, equations, dg, cache, i, j, element)
        (zeta + f)^2 / h
    end
end
