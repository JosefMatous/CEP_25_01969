"""
Calculates the hydrodynamic forces acting on the cable

    F = forces(q, q_dot, mdl)
Calculates the hydrodynamic forces acting on the cable.
Returns a `Vector` of `Pair`s, where each pair represents one hydrodynamic force.

The forces are calculated by solving the following integral:
```
    F = - ∫_0^L (c_u*v_parallel(s) + c_s*v_perpendicular(s)) ds
```
where `v_parallel` and `v_perpendicular` are the components of the relative 
velocity vector that are parallel and perpendicular to the cable's orientation, 
respectively.
The relative velocity vector is given by `v_rel(s) = v(s) - V`, where `v(s)` is
the velocity of the cable at a given point `s`, and `V` is the fluid flow 
velocity (i.e., ocean current).
"""
function EulerLagrangeX.forces(q::Rn, q_dot::Rn, t::Real, mdl::Cable)
    N_q = length(q)

    vector_length = isa(mdl.L, Rn)
    vector_c_u = isa(mdl.c_u, Rn)
    vector_c_s = isa(mdl.c_s, Rn)

    v_begin = q_dot[1:2] - ocean_current(t, mdl.V) # Relative velocity at the beginning of the segment
    cθ = cos(q[3]); sθ = sin(q[3]);
    R = [cθ -sθ; sθ cθ] # rotation matrix
    L = vector_length ? mdl.L[1] : mdl.L
    v_end = v_begin + q_dot[3]*L*[-sθ; cθ] # Relative velocity at the end of the segment

    # Midpoint Jacobian
    T = promote_type(eltype(q), eltype(q_dot), typeof(t))
    J_mid = zeros(T, N_q, 2)
    J_mid[1,1] = 1
    J_mid[2,2] = 1
    J_mid[3,1] = -0.5*L*sθ
    J_mid[3,2] =  0.5*L*cθ

    # Force on the first segment
    ν_begin = R' * v_begin
    ν_end = R' * v_end
    c_u_i = vector_c_u ? mdl.c_u[1] : mdl.c_u
    c_s_i = vector_c_s ? mdl.c_s[1] : mdl.c_s
    # The segment velocities are piecewise linear
    #  => The integral `- ∫_0^L (c_u*v_parallel(s) + c_s*v_perpendicular(s)) ds`
    # can be simplified to: `-L/2 * (c_u*(v_parallel(0)+v_parallel(L)) + c_s*(v_perpendicular(0)+v_perpendicular(L)))`
    τ = -0.5*L * [
        c_u_i*(ν_begin[1] + ν_end[1]);
        c_s_i*(ν_begin[2] + ν_end[2])
        ]
    # Transform the forces back to the inertial frame:
    F = J_mid * (R*τ)

    # Forces on the remaining segments
    for i = 4:N_q
        v_begin = v_end
        J_mid[i-1,1] -= 0.5*L*sθ
        J_mid[i-1,2] += 0.5*L*cθ

        cθ = cos(q[i]); sθ = sin(q[i]);
        R[1,1] = cθ; R[1,2] = -sθ; R[2,1] = sθ; R[2,2] = cθ # rotation matrix
        L = vector_length ? mdl.L[i-2] : mdl.L
        v_end += q_dot[i]*L*[-sθ; cθ]

        J_mid[i,1] = -0.5*L*sθ
        J_mid[i,2] =  0.5*L*cθ

        ν_begin = R' * v_begin
        ν_end = R' * v_end
        c_u_i = vector_c_u ? mdl.c_u[i-2] : mdl.c_u
        c_s_i = vector_c_s ? mdl.c_s[i-2] : mdl.c_s

        τ[1] = -0.5*L * c_u_i*(ν_begin[1] + ν_end[1])
        τ[2] = -0.5*L * c_s_i*(ν_begin[2] + ν_end[2])

        F += J_mid * (R*τ)
    end

    return F
end
