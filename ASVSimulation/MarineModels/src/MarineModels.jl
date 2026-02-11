module MarineModels

    import EulerLagrangeX
    import EulerLagrangeX: Rn, Rmxn
    import DocStringExtensions: FIELDS
    import StaticArrays: @SMatrix, @SVector, @MVector
    import LinearAlgebra

    include("constants.jl")
    include("utility.jl")
    include("model.jl")
    include("types.jl")
    include("so3.jl")    
    include("asv_point.jl")

end # module MarineModels
