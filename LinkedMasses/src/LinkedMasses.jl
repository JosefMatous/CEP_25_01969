module LinkedMasses

    import ForwardDiff
    import DocStringExtensions: FIELDS
    using LinearAlgebra
    using Plots, ColorSchemes
    using Printf
    using LaTeXStrings
    using DifferentialEquations
    using CSV, DataFrames

    const Rn = AbstractVector{<:Real}

    include("model.jl")
    include("simstate.jl")
    include("los.jl")
    include("control.jl")
    include("control_lin.jl")
    include("control_lin_adaptive.jl")
    include("simulation.jl")
    include("visualization.jl")

end # module LinkedMasses
