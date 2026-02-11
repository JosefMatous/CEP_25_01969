export odefun

"""
Ordinary differential equations of the system

    M, b = odefun(q, q_dot, t, elsys)
Calculates the mass matrix `M` and the vector of generalized forces `b` that 
represent the differential equation of the system, i.e.: `M * q_ddot == b`.

    M, b, x_i_dot = odefun(q, q_dot, x_i, e, t, celsys)
Calculates a mass matrix `M`, the vector of generalized forces `b`, and the
derivative of the internal controller state `x_i_dot` for a controlled 
Euler-Lagrange system.
"""
function odefun(q::Rn, q_dot::Rn, t::Real, elsys::EulerLagrangeSystem)
    return _rigid_body_dynamics(q, q_dot, t, elsys)
end

function odefun(q::Rn, q_dot::Rn, x_i::Rn, e::Any, t::Real, celsys::ControlledEulerLagrangeSystem)
    M, b = _rigid_body_dynamics(q, q_dot, t, celsys)
    x_i_dot = add_control_law!(b, q, q_dot, x_i, e, t, celsys)
    return M, b, x_i_dot
end

function _rigid_body_dynamics(q::Rn, q_dot::Rn, t::Real, elsys::Union{EulerLagrangeSystem,ControlledEulerLagrangeSystem})
    # Mass matrix and the Coriolis forces
    M, b = _mass_Coriolis(q, q_dot, t, elsys)
    # Add the generalized forces to the right-hand side of the equation
    add_force!(b, q, q_dot, t, elsys)
    return M, b
end

function _mass_Coriolis(q::Rn, q_dot::Rn, t::Real, elsys::Union{EulerLagrangeSystem,ControlledEulerLagrangeSystem})
    # Euler-Lagrange equation
    #  d/dt(∂ℓ/∂q_dot) - ∂ℓ/∂q = 0
    #  => (∂^2ℓ/∂q_dot^2)*q_ddot + ∂^2ℓ/∂q_dot∂t + (∂^2ℓ/∂q_dot∂q)*q_dot - ∂ℓ/∂q = 0
    #  => M*q_ddot = b, where
    #     M = ∂^2ℓ/∂q_dot^2
    #     b = ∂ℓ/∂q - ∂^2ℓ/∂q_dot∂t - (∂^2ℓ/∂q_dot∂q)*q_dot + Q
    N = length(q)
    x = [q; q_dot; t]
    res = DiffResults.HessianResult(x)
    ForwardDiff.hessian!(res, (x) -> lagrangian(x[1:N], x[N+1:2N], x[end], elsys), x)

    H = DiffResults.hessian(res)
    ∇L = DiffResults.jacobian(res)
    M = H[N+1:2N, N+1:2N]
    b = ∇L[1:N] - H[N+1:2N, end] - H[N+1:2N, 1:N] * q_dot

    return M, b
end

function add_force!(b::Rn, q::Rn, q_dot::Rn, t::Real, elsys::Union{EulerLagrangeSystem,ControlledEulerLagrangeSystem})
    Q = forces(q, q_dot, t, elsys)
    if isa(Q, Rn)
        b .+= Q
    else
        Q_fcn = (q) -> forces(q, q_dot, t, elsys)
        add_positional_force!(b, Q_fcn, q, Q)
    end
end

function add_control_law!(b::Rn, q::Rn, q_dot::Rn, x_i::Rn, e::Any, t::Real, celsys::ControlledEulerLagrangeSystem)
    τ, x_i_dot = control_law(q, q_dot, x_i, e, t, celsys)
    if isa(τ, Rn)
        b .+= τ
    else
        τ_fcn = (q) -> control_law(q, q_dot, x_i, e, t, celsys)[1]
        add_positional_force!(b, τ_fcn, q, τ)
    end
    return x_i_dot
end

function add_positional_force!(b::Rn, Q_fcn::Function, q::Rn, Q)
    Q_position = (q) -> vcat((Q_i.first for Q_i ∈ Q_fcn(q))...)
    J = ForwardDiff.jacobian(Q_position, q)

    counter = 0
    for Q_i ∈ Q
        N_i = length(Q_i.first)
        b .+= J[counter .+ (1:N_i), :]' * Q_i.second
        counter += N_i
    end
end
