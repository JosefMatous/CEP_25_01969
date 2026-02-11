export SimulationState
export get_simulation_state, vectorize

struct SimulationState
    "State of the ASV (position, velocity, internal states)"
    asv_state::Tuple{Rn, Rn, Rn}
    "State of the cable model (position, velocity)"
    cable_state::Tuple{Rn, Rn}
    "State of the payload model (position, velocity)"
    payload_state::Tuple{Rn, Rn}
end

function get_simulation_state(x::Rn, sim::SimulationParams)
    # Sanity checks
    @assert EulerLagrangeX.get_num_coordinates(sim.asv) == 3 "ASV model must have 3 coordinates"
    N_asv = 3
    N_cable = EulerLagrangeX.get_num_coordinates(sim.cable)
    @assert EulerLagrangeX.get_num_coordinates(sim.payload) == 2 "Payload model must have 2 coordinates"
    N_payload = 2
    N_i_asv = EulerLagrangeX.get_num_istates(sim.asv)
    @assert length(x) == 2(N_asv + N_cable + N_payload) + N_i_asv "State vector length does not match the combined model dimensions"

    # Generalized coordinates
    q_asv = x[1:N_asv]
    counter = N_asv
    q_cable = x[(1:N_cable) .+ counter]
    counter += N_cable
    q_payload = x[(1:N_payload) .+ counter]
    counter += N_payload

    # Generalized velocities
    q_dot_asv = x[(1:N_asv) .+ counter]
    counter += N_asv
    q_dot_cable = x[(1:N_cable) .+ counter]
    counter += N_cable
    q_dot_payload = x[(1:N_payload) .+ counter]
    counter += N_payload

    # Internal states
    x_i_asv = x[(1:N_i_asv) .+ counter]

    return SimulationState(
        (q_asv, q_dot_asv, x_i_asv),
        (q_cable, q_dot_cable),
        (q_payload, q_dot_payload)
    )
end

function vectorize(state::SimulationState)
    return vcat(
        state.asv_state[1],
        state.cable_state[1],
        state.payload_state[1],
        state.asv_state[2],
        state.cable_state[2],
        state.payload_state[2],
        state.asv_state[3]
    )
end
