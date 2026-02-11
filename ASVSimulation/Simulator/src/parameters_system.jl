export default_asv_model
export default_cable_model
export default_payload_model
export default_constant_current
export default_varying_current
export default_perturbation

function default_constant_current()
    return [0.5, 0.5]
end

function default_varying_current(t::Real)
    arg = -2π * t / 50.0
    return [0.5 + 0.25 * cos(arg), 0.5 + 0.25 * sin(arg)]
end

function default_perturbation(t::Real, ::Rn)
    arg = 2π * t / 30.0
    return 0.5 * [cos(arg), sin(arg)]
end

function default_asv_model(; V_c=default_constant_current(), μ=zeros(2))
    m = 1000.0
    I_z = 4500.0
    xG = 2.0
    Xu = 500.0
    Yv = 2000.0
    Nr = 1500.0

    Dul = 500.0
    Dvl = 2000.0
    Drl = 4000.0
    Duq = Dul / 2
    Dvq = Dvl / 2
    Drq = Drl / 2

    return MarineModels.ASVModel(
        m, I_z, xG,
        Xu, Yv, Nr, 0.0,
        Dul, Dvl, Drl,
        Duq, Dvq, Drq,
        V_c,
        μ
    )
end

function default_cable_model(; N=8, V_c=default_constant_current(), μ=zeros(2))
    ϱ = 0.5 # mass per unit length
    L = 100.0 / N # length of each segment
    m = ϱ * L # mass of each segment
    c_u = 1.0 # longitudinal drag coefficient
    c_v = 1.0 # transverse drag coefficient

    return CableModels.Cable(N, L, m, c_u, c_v, V_c, μ)        
end

function default_payload_model(; V_c=default_constant_current(), μ=zeros(2))
    m = 250.0
    d = 250.0
    return PayloadModels.Payload(m, d, 0.0, V_c, μ)
end
