export ASVAdaptiveController

"""
    ASVAdaptiveController

Combines line-of-sight guidance with an adaptive linearizing controller for the pendulum dynamics, along with a proportional-derivative heading controller. The internal state includes the LOS path parameter `s` and the adaptive controller parameters `ζ`.

**Fields**
$(FIELDS)
"""
struct ASVAdaptiveController <: MarineModels.Controller
    "Line-of-sight parameters"
    los::LOSParameters
    "Desired path"
    path::Function
    "Adaptive linearizing controller parameters"
    ctrl::AdaptiveLinearizingController
    "Nominal cable length"
    L::Real
    "Heading control proportional gain"
    k_ψ::Real
    "Heading control derivative gain"
    k_r::Real
end

MarineModels.get_num_istates(::ASVAdaptiveController) = parameter_vector_length + 1 # +1 for the LOS path parameter `s`

function MarineModels.control_law(q::Rn, q_dot::Rn, x_i::Rn, e::Any, t::Real, mdl::MarineModels.ASVModel, ctrl::ASVAdaptiveController)
    # Position and velocity of the towed object
    if isa(e, Rn)
        p_tow = e[1:2]
        v_tow = e[3:4]
    elseif isa(e, Tuple) && length(e) == 2
        p_tow, v_tow = e
    else
        error("External factors `e` must be the position and velocity of the towed object, either as a vector of length 4 or as a tuple of two vectors of length 2.")
    end

    # Pendulum state
    p_asv = q[1:2]
    v_asv = q_dot[1:2]
    ψ = q[3]
    r = q_dot[3]
    θ = atanv(p_tow - p_asv)
    θ_dot = ForwardDiff.gradient(atanv, p_tow - p_asv) ⋅ (v_tow - v_asv)
    x_pendulum = [p_asv; θ; v_asv; θ_dot]

    # Internal state
    s = x_i[1]
    ζ = x_i[2:end]

    # Control point
    ε = ctrl.ctrl.ε
    p_ε = (1 - ε) * p_asv + ε * p_tow
    v_ε = (1 - ε) * v_asv + ε * v_tow
    L = norm(p_tow - p_asv) # Effective cable length

    # LOS guidance
    v_LOS, s_dot = line_of_sight(p_ε, s, ctrl.path, ctrl.los)
    v_LOS_fcn = (x_) -> line_of_sight(x_[1:2], x_[3], ctrl.path, ctrl.los)[1]
    ψ_ref = atanv(v_LOS)
    v_LOS_dot = ForwardDiff.jacobian(v_LOS_fcn, [p_ε; s]) * [v_ε; s_dot]
    r_ref = ForwardDiff.gradient(atanv, v_LOS) ⋅ v_LOS_dot

    # Adaptive linearizing control
    f, ζ_dot = adaptive_linearizing_controller(x_pendulum, ζ, v_LOS, v_LOS_dot, L, ctrl.ctrl)

    # Heading control
    ψ_error = mod2pi(ψ - ψ_ref + π) - π
    r_error = r - r_ref
    f_ψ = -ctrl.k_ψ * ψ_error - ctrl.k_r * r_error

    # Suppose the mass of the ASV is known. We can cancel out the Coriolis and added mass effects:
    u = [f; f_ψ]
    M, b = EulerLagrangeX._mass_Coriolis(q, q_dot, t, mdl)
    m = mdl.m
    τ = (M * u) ./ m - b

    # Outputs    
    x_i_dot = [s_dot; ζ_dot]
    return τ, x_i_dot
end
