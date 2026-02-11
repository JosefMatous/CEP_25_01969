module Simulator

    import DocStringExtensions: FIELDS
    import ForwardDiff
    import DifferentialEquations
    import Plots

    import EulerLagrangeX: Rn
    import EulerLagrangeX
    import CableModels
    import PayloadModels
    import MarineModels
    import AdaptiveController

    include("simulation.jl")
    include("simstate.jl")
    include("initial.jl")
    include("parameters_control.jl")
    include("parameters_system.jl")
    include("parameters_simulation.jl")
    include("ode.jl")
    include("get_info.jl")
    include("visualization.jl")

end # module Simulator
