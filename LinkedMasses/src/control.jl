export Controller
export closed_loop_ode!

"Abstract type for controllers of linked masses."
abstract type Controller end

"Returns the linked masses model"
function _get_model(ctrl::Controller)
    try
        return ctrl.model
    catch
        error("Controller does not have a model field.")
    end
end

"Returns the path function"
function _get_path(ctrl::Controller)
    try
        return ctrl.path_fcn
    catch
        error("Controller does not have a path_fcn field.")
    end    
end

"Returns the LOS parameters"
function _get_los(ctrl::Controller)
    try
        return ctrl.los
    catch
        error("Controller does not have a los field.")
    end
end

"Returns the ocean current at the given position and time"
function _get_ocean_current(ctrl::Controller, time::Real)
    try
        return _evaluate_ocean_current(ctrl.V_c, time)
    catch
        error("Controller does not have a V_c field.")
    end
end

_evaluate_ocean_current(V_c::Function, time::Real) = V_c(time)
_evaluate_ocean_current(V_c::Rn, ::Real) = V_c

"""
    closed_loop_ode!(dx, x, params, t)

Closed-loop in-place ODEs for the linked masses system.
"""
closed_loop_ode!(dx::Rn, x::Rn, params::Controller, t::Real) = error("closed_loop_ode! not implemented for $(typeof(params))")
