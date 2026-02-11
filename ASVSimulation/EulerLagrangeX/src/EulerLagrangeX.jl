module EulerLagrangeX

    import DocStringExtensions: FIELDS
    import ForwardDiff, DiffResults

    include("types.jl")
    include("functions.jl")
    include("ode.jl")
    include("object.jl")

end # module EulerLagrange
