export default_initial_state

function default_initial_state(; N_cable=8)
    q_asv = zeros(3)
    q_dot_asv = [1.0, 0.0, 0.0]
    p0_cable = q_asv[1:2]
    θ_cable = range(-π, stop=0.0, length=N_cable)
    q_cable = vcat(p0_cable, collect(θ_cable))
    q_dot_cable = [q_dot_asv[1:2]; zeros(N_cable)]
    L = 100.0 / N_cable
    p_payload = p0_cable + L * [sum(cos.(θ_cable)), sum(sin.(θ_cable))]
    v_payload = q_dot_asv[1:2]

    # Controller internal states
    s0 = 0.0       # path parameter
    # Parameter estimates
    k_err = 2
    m0 = 1000.0
    m = 250.0 * k_err
    c0 = 500.0
    c = 250.0 / k_err
    L = 100.0
    V_c_hat = zeros(2)
    ε = 0.7
    ζ0 = AdaptiveController.parameter_vector(
        AdaptiveController.LinkedMassesParameters(N_cable*L, m0, m, c0, c),
        ε,
        V_c_hat
    )
    x_i_asv = [s0; ζ0]
    
    return SimulationState(
        (q_asv, q_dot_asv, x_i_asv),
        (q_cable, q_dot_cable),
        (p_payload, v_payload)
    )
end
