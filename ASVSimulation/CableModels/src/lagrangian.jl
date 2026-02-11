function EulerLagrangeX.lagrangian(q::Rn, q_dot::Rn, ::Real, mdl::Cable)
    N_q = length(q)

    vector_mass = isa(mdl.m, Rn)
    vector_length = isa(mdl.L, Rn)

    m2 = 0.5 * (vector_mass ? mdl.m[1] : mdl.m) # mass (halved)
    L2 = 0.5 * (vector_length ? mdl.L[1] : mdl.L) # length (halved)
    J2 = m2 * L2^2 / 3 # moment of inertia (halved) (J = 1/12 * m * L^2 ==> J/2 = 1/3 * (0.5*m) * (0.5*L)^2)
    
    cθ = cos(q[3])
    sθ = sin(q[3])
    v_midpoint = q_dot[1:2]
    v_midpoint[1] -= L2*q_dot[3]*sθ
    v_midpoint[2] += L2*q_dot[3]*cθ
    ℓ = m2*LinearAlgebra.dot(v_midpoint, v_midpoint) + J2*q_dot[3]^2 # kinetic energy (T = 1/2*m*|v|^2 + 1/2*J*ω^2)
    for i = 4:N_q
        v_midpoint[1] -= L2*q_dot[i-1]*sθ
        v_midpoint[2] += L2*q_dot[i-1]*cθ

        m2 = 0.5 * (vector_mass ? mdl.m[i-2] : mdl.m) # mass (halved)
        L2 = 0.5 * (vector_length ? mdl.L[i-2] : mdl.L) # length (halved)
        J2 = m2 * L2^2 / 3 # moment of inertia (halved)
        cθ = cos(q[i])
        sθ = sin(q[i])

        v_midpoint[1] -= L2*q_dot[i]*sθ
        v_midpoint[2] += L2*q_dot[i]*cθ

        ℓ += m2*LinearAlgebra.dot(v_midpoint, v_midpoint) + J2*q_dot[i]^2
    end

    return ℓ
end
