export default_line_of_sight
export default_path
export default_adaptive_controller
export default_asv_controller

function default_line_of_sight()
    U = 3.0
    Δ = 15.0
    return AdaptiveController.LOSParameters(U, Δ)
end

function default_path(s::Real)
    R = 100.0
    ϕ0 = -π / 2
    return R * [cos(s + ϕ0), sin(s + ϕ0)] + [20.0, 100.0]
end

function default_adaptive_controller()
    ε = 0.7
    k_v = 0.5
    γ = 5.0
    return AdaptiveController.AdaptiveLinearizingController(k_v, γ, ε)
end

function default_asv_controller()
    los = default_line_of_sight()
    ctrl = default_adaptive_controller()

    k_ψ = 1.0e4
    k_r = 2.0e4

    return AdaptiveController.ASVAdaptiveController(
        los,
        default_path,
        ctrl,
        k_ψ,
        k_r
    )
end
