export Controller, control_law, get_num_istates
export ControlledASV

"Abstract type representing a controller for a marine vehicle"
abstract type Controller end

"""
Performs the control law

    τ, x_i_dot = control_law(q, q_dot, x_i, e, t, mdl, ctrl)
Returns the control forces, `τ`, and the derivative of the internal state, `x_i_dot`.

Arguments:
  - `q`: Generalized coordinates
  - `q_dot`: Generalized velocities
  - `x_i`: Internal state of the controller
  - `e`: External factors (e.g., sensor measurements, reference signals, etc. that cannot be expressed as a function of `q`, `q_dot`, `x_i`, or `t`)
  - `t`: Simulation time
  - `mdl`: ASV model
  - `ctrl`: Controller
"""
function control_law(q::Rn, q_dot::Rn, x_i::Rn, e::Any, t::Real, mdl::ASVModel, ctrl::Controller)
    error("Control law for a vehicle of type $(typeof(mdl)) and a controller of type $(typeof(ctrl)) undefined")
end

"Returns the number of internal states of the controller"
get_num_istates(::Controller) = error("Number of internal states for a controller of type $(typeof(ctrl)) undefined")

"""
A structure encapsulating a marine vehicle model and a controller

**Fields**
$(FIELDS)
"""
struct ControlledASV <: EulerLagrangeX.ControlledEulerLagrangeSystem
    mdl::ASVModel
    ctrl::Controller
end

EulerLagrangeX.get_num_coordinates(sys::ControlledASV) = EulerLagrangeX.get_num_coordinates(sys.mdl)
EulerLagrangeX.get_num_istates(sys::ControlledASV) = get_num_istates(sys.ctrl)
EulerLagrangeX.lagrangian(q::Rn, q_dot::Rn, t::Real, sys::ControlledASV) = EulerLagrangeX.lagrangian(q, q_dot, t, sys.mdl)
EulerLagrangeX.control_law(q::Rn, q_dot::Rn, x_i::Rn, e::Any, t::Real, sys::ControlledASV) = control_law(q, q_dot, x_i, e, t, sys.mdl, sys.ctrl)
EulerLagrangeX.forces(q::Rn, q_dot::Rn, t::Real, sys::ControlledASV) = EulerLagrangeX.forces(q, q_dot, t::Real, sys.mdl)
