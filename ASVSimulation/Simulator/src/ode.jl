export closed_loop_ode

function closed_loop_ode(x::Rn, sim::SimulationParams, t::Real)
    # Get the current simulation state
    state = get_simulation_state(x, sim)

    # System matrices
    N_asv = 3
    N_cable = EulerLagrangeX.get_num_coordinates(sim.cable)
    N_payload = 2
    N_i_asv = EulerLagrangeX.get_num_istates(sim.asv)

    # Euler-Lagrange equations
    T = promote_type(eltype(x), eltype(t))
    # ASV
    M_asv, b_asv, x_i_dot = EulerLagrangeX.odefun(
        state.asv_state[1], state.asv_state[2], state.asv_state[3], state.payload_state,
        t, sim.asv)
    # Cable
    M_cable, b_cable = EulerLagrangeX.odefun(
        state.cable_state[1], state.cable_state[2],
        t, sim.cable)
    # Payload
    M_payload, b_payload = EulerLagrangeX.odefun(
        state.payload_state[1], state.payload_state[2],
        t, sim.payload)

    # Constraints
    k_p = sim.k_f^2
    k_v = 2 * sim.k_f
    # 1. ASV-cable attachment point
    asv_obj = EulerLagrangeX.Object(sim.asv, 1)
    cable_obj = EulerLagrangeX.Object(sim.cable, 2)
    asv_point = MarineModels.ASVPoint(asv_obj, [0.0, 0.0])
    cable_begin = CableModels.cable_begin(cable_obj)

    p_asv, v_asv = EulerLagrangeX.get_point(
        state.asv_state[1], state.asv_state[2],
        asv_point)
    J_asv, a0_asv = EulerLagrangeX.get_point_acceleration(
        state.asv_state[1], state.asv_state[2],
        asv_point
    )

    p_cable_begin, v_cable_begin = EulerLagrangeX.get_point(
        state.cable_state[1], state.cable_state[2],
        cable_begin)
    J_cable_begin, a0_cable_begin = EulerLagrangeX.get_point_acceleration(
        state.cable_state[1], state.cable_state[2],
        cable_begin)

    # 2. Cable-payload attachment point
    cable_end = CableModels.cable_end(cable_obj)
    payload_obj = EulerLagrangeX.Object(sim.payload, 3)
    payload_point = PayloadModels.PayloadPoint(payload_obj)

    p_cable_end, v_cable_end = EulerLagrangeX.get_point(
        state.cable_state[1], state.cable_state[2],
        cable_end)
    J_cable_end, a0_cable_end = EulerLagrangeX.get_point_acceleration(
        state.cable_state[1], state.cable_state[2],
        cable_end)

    p_payload, v_payload = EulerLagrangeX.get_point(
        state.payload_state[1], state.payload_state[2],
        payload_point)
    J_payload, a0_payload = EulerLagrangeX.get_point_acceleration(
        state.payload_state[1], state.payload_state[2],
        payload_point)

    # Assemble the full system
    #  M * [q_ddot_asv; q_ddot_cable; q_ddot_payload; λ_asv_cable; λ_cable_payload] == b
    # Mass matrix:
    M = [
                  M_asv               zeros(T, N_asv, N_cable)    zeros(T, N_asv, N_payload)         J_asv'          zeros(T, N_asv, 2);
         zeros(T, N_cable, N_asv)             M_cable            zeros(T, N_cable, N_payload)    -J_cable_begin'        J_cable_end';
        zeros(T, N_payload, N_asv) zeros(T, N_payload, N_cable)            M_payload          zeros(T, N_payload, 2)    -J_payload';
                  J_asv                    -J_cable_begin          zeros(T, 2, N_payload)        zeros(T, 2, 2)        zeros(T, 2, 2);
              zeros(T, 2, N_asv)             J_cable_end                  -J_payload             zeros(T, 2, 2)        zeros(T, 2, 2)        
    ]
    # Right-hand side:
    b = vcat(
        b_asv,
        Vector(b_cable), # Convert to Vector, otherwise the result is a sparse vector
        b_payload,
        -k_p * (p_asv - p_cable_begin) - k_v * (v_asv - v_cable_begin) - (a0_asv - a0_cable_begin),
        -k_p * (p_cable_end - p_payload) - k_v * (v_cable_end - v_payload) - (a0_cable_end - a0_payload)
    )
    # Solve for accelerations and Lagrange multipliers
    z = M \ b

    # Populate the state derivative
    #  dx = [q_dot_asv; q_dot_cable; q_dot_payload; q_ddot_asv; q_ddot_cable; q_ddot_payload; x_i_dot_asv]
    dx = [
        state.asv_state[2]; # q_dot_asv
        state.cable_state[2]; # q_dot_cable
        state.payload_state[2]; # q_dot_payload
        z[1:N_asv]; # q_ddot_asv
        z[(N_asv+1):(N_asv+N_cable)]; # q_ddot_cable
        z[(N_asv+N_cable+1):(N_asv+N_cable+N_payload)]; # q_ddot_payload
        x_i_dot; # x_i_dot_asv
    ]

    return dx
end
