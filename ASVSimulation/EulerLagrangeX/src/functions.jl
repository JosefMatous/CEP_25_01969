export get_num_coordinates, get_num_istates
export lagrangian, forces, control_law

"""
Returns the number of generalized coordinates of the Euler-Lagrange system
"""
get_num_coordinates(sys::Union{EulerLagrangeSystem, ControlledEulerLagrangeSystem}) = 0

"""
Returns the number of internal states of a controlled Euler-Lagrange system
"""
get_num_istates(sys::Union{EulerLagrangeSystem, ControlledEulerLagrangeSystem}) = 0

"""
Calculates the Lagrangian of the system

    ℓ = lagrangian(q, q_dot, t, sys)
Returns the Lagrangian of the system for the given configuration. The Lagrangian is given by
`ℓ = T(q,q_dot) - P(q)`, where `T` and `P` is the kinetic and potential energy
of the system, respectively.

Arguments:
  - `q`: Generalized coordinates
  - `q_dot`: Generalized velocities
  - 't': Simulation time
  - `sys`: System 
"""
function lagrangian(q::Rn, q_dot::Rn, t::Real, sys::Union{EulerLagrangeSystem, ControlledEulerLagrangeSystem})
    error("Missing the implementation of lagrangian for $(typeof(sys))")
end

"""
Returns the forces acting on the system

    F = forces(q, q_dot, t, sys)
Calculates the forces acting on the system (`F`). `F` can either be
  1. A vector of generalized forces (i.e., a vector of real values of the same length as `q`)
  2. A list (i.e., a `Vector` or a `Tuple`) of `Pair{Vector, Vector}`.
Each `Pair` represents a position where the force is applied and the magnitude of the force,
i.e., `F[i] = p_i => f_i`, where `p_i` is the position (N-D vector), and `f_i` is the force (N-D vector),
where `N` is the number of dimensions.

Arguments:
  - `q`: Generalized coordinates
  - `q_dot`: Generalized velocities
  - `t`: Simulation time
  - `sys`: System 
"""
function forces(q::Rn, q_dot::Rn, t::Real, sys::Union{EulerLagrangeSystem, ControlledEulerLagrangeSystem})
    error("Missing the implementation of forces for $(typeof(sys))")
end

"""
Performs the control law

    τ, x_i_dot = control_law(q, q_dot, x_i, e, t, sys)
Returns the control forces, `τ`, and the derivative of the internal state, `x_i_dot`.

`τ` can either be
  1. A vector of generalized forces (i.e., a vector of real values of the same length as `q`)
  2. A list (i.e., a `Vector` or a `Tuple`) of `Pair{Vector, Vector}`.
Each `Pair` represents a position where the force is applied and the magnitude of the force,
i.e., `τ[i] = p_i => f_i`, where `p_i` is the position (N-D vector), and `f_i` is the force (N-D vector),
where `N` is the number of dimensions.

`x_i_dot` is a vector of real values of the same size as `x_i`.

Arguments:
  - `q`: Generalized coordinates
  - `q_dot`: Generalized velocities
  - `x_i`: Internal state of the controller
  - `e`: External factors (e.g., sensor measurements, reference signals, etc. that cannot be expressed as a function of `q`, `q_dot`, `x_i`, or `t`)
  - `t`: Simulation time
  - `sys`: System 

See also: `forces`
"""
function control_law(q::Rn, q_dot::Rn, x_i::Rn, e::Any, t::Real, sys::ControlledEulerLagrangeSystem)
    error("Missing the implementation of control_law for $(typeof(sys))")
end

function control_law(q::Rn, q_dot::Rn, x_i::Rn, e::Any, t::Real, sys::EulerLagrangeSystem)
    return [], []
end
