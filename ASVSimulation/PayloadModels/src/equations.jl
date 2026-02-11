function EulerLagrangeX.lagrangian(
    ::Rn, 
    q_dot::Rn, 
    ::Real, 
    sys::Payload
) 
    m = sys.m

    # Kinetic energy
    T = 0.5 * m * sum((q_dot).^2)
    return T 
end

function EulerLagrangeX.forces(
    q::Rn, 
    q_dot::Rn, 
    t::Real, 
    mdl::Payload
)
    
    V_c = ocean_current(t, mdl.V)
    μ = perturbation(t, q, mdl.μ)

    v_rel = q_dot - V_c
    U = sqrt(sum(v_rel.^2))

    F_lin = -mdl.d_l * v_rel
    F_quad = -mdl.d_q * U * v_rel
    Q_pert = mdl.m * μ

    Q = F_lin + F_quad + Q_pert
    return Q
end
