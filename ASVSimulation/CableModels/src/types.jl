export Cable

"""
Parameters of a rigid-body cable model

A planar (2D) model of a cable consisting of multiple rigid-body segments. 
The model has `N+2` degrees of freedom, where `N` is the number of segments.
The generalized coordinates of the system are: `q = [p0; θ]`, where `p0 ∈ R^2`
is the position of the first segment, and `θ ∈ R^N` is the orientation of the
segment (relative to the inertial coordinate frame).

Fields:
$(FIELDS)

Constructors:
    `Cable(N, L, m, c_u, c_s, V=[0,0], μ=[0,0])`
Creates a cable with `N` segments. `L`, `m`, `c_u`, and `c_s` can be vectors or scalars.
If these parameters are scalar, then it is assumed that the value of these parameters
is identical for all segments.

    `Cable(L, m, c_u, c_s, V=[0,0], μ=[0,0])`
Deduces the number of segments from `L`, i.e., `N = length(L)`.
"""
struct Cable <: EulerLagrangeX.EulerLagrangeSystem
    "Number of cable segments"
    N::Integer
    "Length of cable segments"
    L::Union{Real, Rn}
    "Mass of cable segments"
    m::Union{Real, Rn}
    "Longitudinal damping"
    c_u::Union{Real, Rn}
    "Lateral damping"
    c_s::Union{Real, Rn}
    "Fluid flow velocity (ocean current) (a 2D vector or a function that accepts one scalar - time - and returns a 2D vector)"
    V::Union{Rn, Function}
    "Perturbation (a 2D vector or a function that accepts a scalar and a 2D vector - time and position - and returns a 2D vector)"
    μ::Union{Rn, Function}

    function Cable(N::Integer, L::Union{Real, Rn}, m::Union{Real, Rn}, c_u::Union{Real, Rn}, c_s::Union{Real, Rn}, V::Union{Rn, Function}=[0,0], μ::Union{Rn, Function}=[0,0])
        if isa(L, Rn)
            @assert (length(L) == N) "Length of L must be equal to N"
        end
        if isa(m, Rn)
            @assert (length(m) == N) "Length of m must be equal to N"
        end
        if isa(c_u, Rn)
            @assert (length(c_u) == N) "Length of c_u must be equal to N"
        end
        if isa(c_s, Rn)
            @assert (length(c_s) == N) "Length of c_s must be equal to N"
        end
        if isa(V, Rn)
            @assert (length(V) == 2) "V must be a 2D vector"
        else
            try
                V0 = V(0.0)
                @assert isa(V0, Rn) && length(V0) == 2
            catch
                error("V must be a function that accepts one scalar argument - time - and returns a 2D vector")
            end
        end
        if isa(μ, Rn)
            @assert (length(μ) == 2) "μ must be a 2D vector"
        else
            try
                μ0 = μ(0.0, [0.0, 0.0])
                @assert isa(μ0, Rn) && length(μ0) == 2
            catch
                error("μ must be a function that accepts two arguments - time and position - and returns a 2D vector")
            end
        end

        return new(N, L, m, c_u, c_s, V, μ)
    end
end

function Cable(L::Rn, m::Union{Real, Rn}, c_u::Union{Real, Rn}, c_s::Union{Real, Rn}, V::Union{Rn,Function}=[0;0], μ::Union{Rn, Function}=[0,0])
    N = length(L)
    return Cable(N, L, m, c_u, c_s, V, μ)
end

EulerLagrangeX.get_num_coordinates(mdl::Cable) = 2 + mdl.N

ocean_current(::Real, V::Rn) = V
ocean_current(t::Real, V::Function) = V(t)

perturbation(::Real, ::Rn, μ::Rn) = μ
perturbation(t::Real, x::Rn, μ::Function) = μ(t, x)
