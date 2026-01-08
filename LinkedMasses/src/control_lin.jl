export LinearizingController

"""
    LinearizingController

Line-of-sight + feedback-linearizing controller for two linked masses.

**Fields**
$(FIELDS)
"""
struct LinearizingController <: Controller
    "Parameters of the linked masses system"
    model::LinkedMassesParameters
    "Line-of-sight guidance parameters"
    los::LOSParameters
    "Path function"
    path_fcn::Function
    "Velocity feedback gain"
    k_v::Real
    "Control point coefficient"
    ε::Real
    "Ocean current"
    V_c::Union{Function, Rn}
end

function closed_loop_ode!(dx::Rn, x::Rn, params::LinearizingController, t::Real)
    # Unpack the state
    p0 = x[1:2]
    θ = x[3]
    v0 = x[4:5]
    θ_dot = x[6]
    s = x[7]
    # Unpack the parameters
    lm = params.model
    los = params.los
    path_fcn = params.path_fcn
    k_v = params.k_v
    ε = params.ε
    V_c = _evaluate_ocean_current(params.V_c, t)

    # Control point
    Γ = [cos(θ), sin(θ)]
    dΓ = [-sin(θ), cos(θ)]
    p = p0 + ε * lm.L * Γ
    v = v0 + ε * lm.L * θ_dot * dΓ

    # LOS guidance
    v_LOS, s_dot = line_of_sight(p, s, path_fcn, los)
    v_LOS_dot = ForwardDiff.jacobian(
        (ξ) -> line_of_sight(ξ[1:2], ξ[3], path_fcn, los)[1],
        [p; s]
    ) * [v; s_dot]

    # Control law
    f = linearizing_controller(x, v_LOS, v_LOS_dot, lm, ε, k_v, V_c)

    # System equations
    dx[1:6] = linked_masses_ode(x[1:6], f, lm, V_c)
    dx[7] = s_dot
    nothing
end

"""
    linearizing_controller(x, v_ref, v_ref_dot, params, ε, k_v, V_c=zeros(2))

Output-linearizing velocity feedback controller for two linked masses.

**Arguments**
- `x::Rn`: State vector
- `v_ref::Rn`: Reference velocity (2D vector).
- `v_ref_dot::Rn`: Time derivative of the reference velocity (2D vector).
- `params::LinkedMassesParameters`: Structure containing system parameters.
- `ε::Real`: Control point coefficient.
- `k_v::Real`: Velocity feedback gain.
- `V_c::Rn=zeros(2)`: Fluid flow velocity (default is zero).
"""
function linearizing_controller(x::Rn, v_ref::Rn, v_ref_dot::Rn, params::LinkedMassesParameters, ε::Real, k_v::Real, V_c::Rn=zeros(2))
    # Unpack parameters
    L = params.L
    m0 = params.m0
    m = params.m
    c0 = params.c0
    c = params.c

    # Unpack states
    θ = x[3]
    v0 = x[4:5]
    θ_dot = x[6]

    Γ = [cos(θ), sin(θ)]
    dΓ = [-sin(θ), cos(θ)]

    # Virtual output
    #p = p0 + ε * L * Γ
    v = v0 + ε * L * θ_dot * dΓ
    v_err = v - v_ref

    # Hydrodynamic forces
    F0 = -c0 * (v0 - V_c) # Damping force on the first mass
    v1 = v0 + L * θ_dot * dΓ # Velocity of the second mass
    F = -c * (v1 - V_c) # Damping force on the second mass

    # Output derivative
    #   M * v_dot = Q + f
    J = dΓ * dΓ' # projection matrix
    M = (m0 + m) * I(2) - (m + ε*m0/(ε - 1)) * J # Control point mass matrix
    Q = (m*(1 - ε) - m0*ε) * L * θ_dot^2 * Γ +
            (I(2) - (1 + m0*ε/(m*(ε - 1))) * J) * F +
            F0 # Control point forces
    M_dot = θ_dot * (m + ε*m0/(ε - 1)) * (Γ * dΓ' + dΓ * Γ') # Control point mass matrix derivative
    
    # Lyapunov-based control law
    # V = 1/2 * (v_err' * M * v_err)
    # V_dot = 1/2 * (v_err' * M_dot * v_err) + v_err' * (Q + f) - v_err' * M * v_ref_dot
    f = M * v_ref_dot  -  Q  -  k_v * M * v_err  -  0.5 * (M_dot * v_err)
    return f
end
