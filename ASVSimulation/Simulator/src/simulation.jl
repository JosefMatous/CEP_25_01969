export SimulationParams

struct SimulationParams
    "ASV parameters and controller"
    asv::MarineModels.ControlledASV
    "Cable"
    cable::CableModels.Cable
    "Payload"
    payload::PayloadModels.Payload
    "Constraint enforcement gain"
    k_f::Real
end
