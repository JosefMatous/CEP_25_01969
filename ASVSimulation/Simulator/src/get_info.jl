export ClosedLoopState, get_closed_loop_state
export PlottingState, get_plotting_state

struct ClosedLoopState
    "Path-following error x"
    e_x::Real
    "Path-following error y"
    e_y::Real
    "Velocity error x"
    v_x::Real
    "Velocity error y"
    v_y::Real
    "Cable angular velocities"
    θ_dot::Rn
end

atanv(v::Rn) = atan(v[2], v[1])

function get_closed_loop_state(x::Rn, sim::SimulationParams)
    state = get_simulation_state(x, sim)
    path = sim.asv.ctrl.path
    los_params = sim.asv.ctrl.los
    ε = sim.asv.ctrl.ctrl.ε

    p0 = state.asv_state[1][1:2]
    p = state.payload_state[1]
    v0 = state.asv_state[2][1:2]
    v = state.payload_state[2]

    s = state.asv_state[3][1]

    p_ε = (1 - ε) * p0 + ε * p
    v_ε = (1 - ε) * v0 + ε * v

    p_path = path(s)
    ∂p_path = ForwardDiff.derivative(path, s)
    θ_path = atanv(∂p_path)
    R_path = [cos(θ_path) -sin(θ_path); sin(θ_path) cos(θ_path)]
    e = R_path' * (p_ε - p_path)

    v_ref, _ = AdaptiveController.line_of_sight(p_ε, s, path, los_params)
    v_err = v_ε - v_ref

    θ_dot = state.cable_state[2][3:end]

    return ClosedLoopState(
        e[1],
        e[2],
        v_err[1],
        v_err[2],
        θ_dot
    )
end

struct PlottingState
    "ASV coordinates"
    q_asv::Rn
    "Cable knot x-coordinates"
    x_cable::Rn
    "Cable knot y-coordinates"
    y_cable::Rn
    "Payload coordinates"
    q_payload::Rn
    "Control point"
    p_ε::Rn
    "Path reference"
    p_path::Rn
end

function get_plotting_state(x::Rn, sim::SimulationParams)
    state = get_simulation_state(x, sim)
    path = sim.asv.ctrl.path
    ε = sim.asv.ctrl.ctrl.ε

    q_asv = state.asv_state[1]
    q_cable = state.cable_state[1]
    q_payload = state.payload_state[1]

    p0 = q_asv[1:2]
    p = q_payload[1:2]

    p_ε = (1 - ε) * p0 + ε * p

    s = state.asv_state[3][1]
    p_path = path(s)

    p0_cable = q_cable[1:2]
    θ_cable = q_cable[3:end]
    x_cable, y_cable = CableModels.config2knots(p0_cable, θ_cable, sim.cable.L)

    return PlottingState(
        q_asv,
        x_cable,
        y_cable,
        q_payload,
        p_ε,
        p_path
    )
end
