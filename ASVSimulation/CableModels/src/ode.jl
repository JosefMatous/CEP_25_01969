function EulerLagrangeX.odefun(q::Rn, q_dot::Rn, t::Real, mdl::Cable)
    N_q = length(q)

    vector_mass = isa(mdl.m, Rn)
    vector_length = isa(mdl.L, Rn)

    m = vector_mass ? mdl.m[1] : mdl.m
    L = vector_length ? mdl.L[1] : mdl.L
    J = (m * L^2) / 12 # moment of inertia

    cθ = cos(q[3])
    sθ = sin(q[3])
    ω = q_dot[3]

    # Midpoint velocity
    v_midpoint = q_dot[1:2]
    v_midpoint[1] -= 0.5*L*ω*sθ
    v_midpoint[2] += 0.5*L*ω*cθ

    T = promote_type(eltype(q), eltype(q_dot), typeof(t))

    # Jacobian of midpoint velocity with respect to q
    Jac_q = SparseArrays.spzeros(T, 2, N_q)
    Jac_q[1,3] = -0.5*L*ω*cθ
    Jac_q[2,3] = -0.5*L*ω*sθ

    # Jacobian of midpoint velocity with respect to q_dot
    Jac_q_dot = SparseArrays.spzeros(T, 2, N_q)
    Jac_q_dot[1,1] = 1
    Jac_q_dot[2,2] = 1
    Jac_q_dot[1,3] = -0.5*L*sθ
    Jac_q_dot[2,3] =  0.5*L*cθ

    # Midpoint velocity times the Hessian of midpoint velocity with respect to q and q_dot 
    #  `Hess_1 = ∂^2v_mid[1]/∂q∂q_dot`
    #  `Hess_2 = ∂^2v_mid[2]/∂q∂q_dot`
    # Note that the result is a diagonal matrix, so we only need to store the diagonal entries:
    Hess_1 = SparseArrays.spzeros(T, N_q)
    Hess_1[3] = -0.5*L*cθ
    Hess_2 = SparseArrays.spzeros(T, N_q)
    Hess_2[3] = -0.5*L*sθ

    # Mass matrix
    #  `M = ∂^2ℓ/∂q_dot^2 = ∑_i m_i (∂v_i/∂q_dot)' * ∂v_i/∂q_dot   +   diag(J)`
    M = m * Jac_q_dot' * Jac_q_dot
    M[3,3] += J

    # Coriolis matrix
    #  `C = ∂^2ℓ/∂q∂q_dot = ∑_i m_i ( (∂v_i/∂q_dot)' * ∂v_i/∂q  +  v_i * ∂^2v_i/∂q∂q_dot`
    # C*q_dot
    b = -m * (Jac_q_dot'*(Jac_q*q_dot) + (v_midpoint[1]*Hess_1 + v_midpoint[2]*Hess_2).*q_dot)

    # Potential effects:
    #  `∇_q ℓ = ∑_i m_i (∂v_i/∂q)' * v_i`
    b .+= m * (Jac_q' * v_midpoint)

    # Perturbation
    #  `Q = ∑_i m_i (∂p_i/∂q)' * μ_i`
    p_midpoint = q[1:2]
    p_midpoint[1] += 0.5*L*cθ
    p_midpoint[2] += 0.5*L*sθ
    b .+= m * (Jac_q_dot' * perturbation(t, p_midpoint, mdl.μ))

    # Iterate over the remaining midpoints
    for i = 4:N_q
        v_midpoint[1] -= 0.5*L*ω*sθ
        v_midpoint[2] += 0.5*L*ω*cθ

        p_midpoint[1] += 0.5*L*cθ
        p_midpoint[2] += 0.5*L*sθ

        cθ = cos(q[i])
        sθ = sin(q[i])
        ω = q_dot[i]

        m = vector_mass ? mdl.m[i-2] : mdl.m
        L = vector_length ? mdl.L[i-2] : mdl.L
        J = (m * L^2) / 12 # moment of inertia

        v_midpoint[1] -= 0.5*L*ω*sθ
        v_midpoint[2] += 0.5*L*ω*cθ

        p_midpoint[1] += 0.5*L*cθ
        p_midpoint[2] += 0.5*L*sθ

        Jac_q[:,i-1] .*= 2
        Jac_q[1,i] = -0.5*L*ω*cθ
        Jac_q[2,i] = -0.5*L*ω*sθ

        Jac_q_dot[:,i-1] .*= 2
        Jac_q_dot[1,i] = -0.5*L*sθ
        Jac_q_dot[2,i] =  0.5*L*cθ

        Hess_1[i-1] *= 2
        Hess_1[i] = -0.5*L*cθ
        Hess_2[i-1] *= 2
        Hess_2[i] = -0.5*L*sθ

        M .+= m * Jac_q_dot' * Jac_q_dot
        M[i,i] += J

        b .-= m * (Jac_q_dot'*(Jac_q*q_dot) + (v_midpoint[1]*Hess_1 + v_midpoint[2]*Hess_2).*q_dot)
        b .+= m * (Jac_q' * v_midpoint)
        b .+= m * (Jac_q_dot' * perturbation(t, p_midpoint, mdl.μ))
    end

    # Hydrodynamic forces
    b .+= EulerLagrangeX.forces(q, q_dot, t, mdl)

    return M, b
end
