export line_of_sight, LOSParameters

"""
    struct LOSParameters

A structure to represent the parameters of line-of-sight (LOS) guidance.

**Fields**
- `U::Real`: Speed.
- `Δ::Real`: Lookahead distance.
"""
struct LOSParameters
    "Speed"
    U::Real
    "Lookahead distance"
    Δ::Real
end

atanv(x::Rn) = atan(x[2], x[1])

"""
    line_of_sight(p::Rn, s::Real, path_fcn::Function, params::LOSParameters)

Line-of-sight (LOS) guidance law.

**Arguments**
- `p::Rn`: Position (2D).
- `s::Real`: Path parameter.
- `path_fcn::Function`: Path parametrization.
- `params::LOSParameters`: Line-of-sight parameters.

**Returns**
LOS velocity and path parameter update.
"""
function line_of_sight(p::Rn, s::Real, path_fcn::Function, params::LOSParameters)
    # Unpack parameters
    U = params.U
    Δ = params.Δ

    # Compute the path point and its derivative
    p_path = path_fcn(s)
    p_path_dot = ForwardDiff.derivative(path_fcn, s)
    θ_path = atanv(p_path_dot)
    R_path = [cos(θ_path) -sin(θ_path); sin(θ_path) cos(θ_path)]
    ∇_norm = norm(p_path_dot)

    # Path-following error
    e = R_path' * (p - p_path)
    e_x = e[1]
    e_y = e[2]
    D = sqrt(Δ^2 + e_y^2)

    # LOS guidance law
    v_LOS = U / D * R_path * [Δ, -e_y]

    # Path parameter update
    s_dot = U / ∇_norm * (Δ/D + e_x/sqrt(Δ^2 + e_x^2))

    return v_LOS, s_dot
end
