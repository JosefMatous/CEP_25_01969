
"""
Total derivative of four-quadrant inverse tan.

    dr = atan_dot(y, x, y_dot, x_dot)
Returns the total derivative of `atan(y,x)`, assuming that `dy/dt = y_dot`,
and `dx/dt = x_dot`.
"""
function atan_dot(y::Real, x::Real, y_dot::Real, x_dot::Real)
    return (x*y_dot - y*x_dot) / (x^2 + y^2)
end

"""
Cross product of a 2D vector.
"""
function cross_2D(a::Rn, b::Rn)
    return a[1]*b[2] - a[2]*b[1]
end

ocean_current(::Real, V::Rn) = V
ocean_current(t::Real, V::Function) = V(t)

perturbation(::Real, ::Rn, μ::Rn) = μ
perturbation(t::Real, p::Rn, μ::Function) = μ(t, p)
