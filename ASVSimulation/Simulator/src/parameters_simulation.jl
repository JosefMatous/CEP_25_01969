export default_simulation

function default_simulation(; N=8, V_c=default_constant_current(), μ=zeros(2), k_f=1.0)
    asv_mdl = default_asv_model(V_c=V_c, μ=μ)
    cable = default_cable_model(N=N, V_c=V_c, μ=μ)
    payload = default_payload_model(V_c=V_c, μ=μ)
    ctrl = default_asv_controller()

    asv = MarineModels.ControlledASV(asv_mdl, ctrl)

    return SimulationParams(asv, cable, payload, k_f)
end
