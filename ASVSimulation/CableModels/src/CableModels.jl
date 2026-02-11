module CableModels

    import LinearAlgebra
    import EulerLagrangeX
    import EulerLagrangeX: Rn, Rmxn
    import DocStringExtensions: FIELDS
    import SparseArrays

    include("types.jl")
    include("coordinates.jl")
    include("hydrodynamics.jl")
    include("lagrangian.jl")
    include("ode.jl")
    include("cablepoint.jl")    

end # module CableModel
