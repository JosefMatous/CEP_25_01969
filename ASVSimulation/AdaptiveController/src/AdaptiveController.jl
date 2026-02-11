module AdaptiveController

    import DocStringExtensions: FIELDS
    import ForwardDiff
    import LinearAlgebra: I, norm, ⋅
    import EulerLagrangeX
    import EulerLagrangeX: Rn, Rmxn
    import MarineModels

    include("linked_masses.jl")
    include("los.jl")
    include("control_lin_adaptive.jl")
    include("asv_control.jl")

end # module AdaptiveController
