export AdaptiveLinearizingController, parameter_vector

"""
    AdaptiveLinearizingController

Line-of-sight controller for two linked masses with model reference adaptive control (MRAC).

**Fields**
$(FIELDS)
"""
struct AdaptiveLinearizingController <: Controller
    "Parameters of the linked masses system"
    model::LinkedMassesParameters
    "Line-of-sight guidance parameters"
    los::LOSParameters
    "Path function"
    path_fcn::Function
    "Velocity feedback gain"
    k_v::Real
    "Adaptation gain"
    k_a::Real
    "Control point coefficient"
    ε::Real
    "Ocean current"
    V_c::Union{Function, Rn}
end

function closed_loop_ode!(dx::Rn, x::Rn, params::AdaptiveLinearizingController, t::Real)
    # Unpack the state
    p0 = x[1:2]
    θ = x[3]
    v0 = x[4:5]
    θ_dot = x[6]
    s = x[7]
    ζ_hat = x[8:end]
    # Unpack the parameters
    lm = params.model
    los = params.los
    path_fcn = params.path_fcn
    k_v = params.k_v
    k_a = params.k_a
    ε = params.ε
    V_c = _evaluate_ocean_current(params.V_c, t) # Remark: only used by the ODE function, not the controller

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
    f, ζ_dot = adaptive_linearizing_controller(x, ζ_hat, v_LOS, v_LOS_dot, lm, ε, k_v, k_a)

    dx[1:6] = linked_masses_ode(x[1:6], f, lm, V_c)
    dx[7] = s_dot
    dx[8:end] = ζ_dot
    nothing
end

"""
    adaptive_linearizing_controller(
        x::Rn, 
        v_ref::Rn, 
        v_ref_dot::Rn, 
        params::LinkedMassesParameters, 
        ε::Real, 
        k_v::Real, 
        k_a::Real
    ) -> (f, ζ_dot)

Adaptive controller for a system of linked masses.

**Arguments**
- `x::Rn`: State vector
- `ζ::Rn`: Parameter estimates vector
- `v_ref::Rn`: Reference velocity vector.
- `v_ref_dot::Rn`: Reference acceleration vector.
- `params::LinkedMassesParameters`: Structure containing system parameters (e.g., link length `L`).
- `ε::Real`: Control point coefficient.
- `k_v::Real`: Velocity feedback gain.
- `k_a::Real`: Adaptation gain.

**Returns**
- `f`: Control input vector for the system.
- `ζ_dot`: Time derivative of the parameter estimates (adaptation law).

**Description**
This function implements an adaptive control law for a planar system of two linked masses. The controller uses a regressor matrix and parameter estimates to compute the control input, and updates the parameter estimates using an adaptation law. The control law is designed to track a reference velocity while compensating for unknown system parameters.
For the structure of the parameter vector `ζ`, see `parameter_vector`.
"""
function adaptive_linearizing_controller(x::Rn, ζ::Rn, v_ref::Rn, v_ref_dot::Rn, params::LinkedMassesParameters, ε::Real, k_v::Real, k_a::Real)
    # Unpack parameters
    L = params.L
    # The parameters below are assumed unknown
    #m0 = params.m0
    #m = params.m
    #c0 = params.c0
    #c = params.c

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

    v1 = v0 + L * θ_dot * dΓ # Velocity of the second mass
    J = dΓ * dΓ' # Projection matrix

    # Regressor
    Y = hcat(
        L * θ_dot^2 * Γ   -   θ_dot / (2 * (ε - 1)) * (Γ * dΓ' + dΓ * Γ') * v_err  -  J * (v_ref_dot - k_v * v_err) ./ (ε - 1),
        v0,
        L * θ_dot * dΓ,
        I(2),
        J * v1,
        J,
        k_v * v_err - v_ref_dot
    )

    # Adaptive control law
    f = -Y * ζ
    # Adaptation law
    ζ_dot = k_a * (Y' * v_err)

    return f, ζ_dot
end

"""
    parameter_vector(param::LinkedMassesParameters, ε::Real, V_c::Rn=zeros(2))

Compute the parameter vector `ζ ∈ R^9` for a system of linked masses. 
The parameter vector is chosen such that the system dynamics are affine in `ζ`.
The vector is defined as:
```
ζ = [
    m * (1 - ε) - m0 * ε;
    -(c + c0);
    -c;
    (c + c0) * V_c;
    c * (1 + m0 * ε / (m * (ε - 1)));
    -c * (1 + m0 * ε / (m * (ε - 1))) * V_c;
    m + m0
]
```

**Arguments**
- `param::LinkedMassesParameters`: Struct containing the physical parameters of the linked masses system, including `m0`, `m`, `c0`, and `c`.
- `ε::Real`: A real-valued parameter, typically representing a scaling or perturbation factor in the system.
- `V_c::Rn`: Ocean current.
"""
function parameter_vector(param::LinkedMassesParameters, ε::Real, V_c::Rn=zeros(2))
    # Unpack parameters
    #L = param.L # The only known parameter
    m0 = param.m0
    m = param.m
    c0 = param.c0
    c = param.c
    
    ζ = [
        m * (1 - ε) - m0 * ε;
        -(c + c0);
        -c;
        (c + c0) * V_c;
        c * (1 + m0 * ε / (m * (ε - 1)));
        -c * (1 + m0 * ε / (m * (ε - 1))) * V_c;
        m + m0
    ]
    return ζ
end
