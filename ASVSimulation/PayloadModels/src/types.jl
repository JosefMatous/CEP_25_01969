export Payload, PayloadPoint

struct Payload <: EulerLagrangeX.EulerLagrangeSystem
    "Mass of the payload"
    m::Real
    "Linear damping coefficient"
    d_l::Real
    "Quadratic damping coefficient"
    d_q::Real
    "Fluid flow velocity (ocean current) (a 2D vector or a function that accepts one scalar - time - and returns a 2D vector)"
    V::Union{Rn, Function}
    "Perturbation (a 2D vector or a function that accepts a scalar and a 2D vector - time and position - and returns a 2D vector)"
    μ::Union{Rn, Function}
end

EulerLagrangeX.get_num_coordinates(::Payload) = 2

struct PayloadPoint <: EulerLagrangeX.ObjectPoint
    "The payload"
    obj::EulerLagrangeX.Object{Payload}
end

EulerLagrangeX.get_num_dimensions(::PayloadPoint) = 2

function EulerLagrangeX.get_point(q::Rn, q_dot::Rn, ::PayloadPoint)
    p = copy(q)
    v = copy(q_dot)
    return p, v
end

ocean_current(::Real, V::Rn) = V
ocean_current(t::Real, V::Function) = V(t)

perturbation(::Real, ::Rn, μ::Rn) = μ
perturbation(t::Real, x::Rn, μ::Function) = μ(t, x)
