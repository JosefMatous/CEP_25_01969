export get_endpoint, get_midpoints
export config2knots, knots2config

"""
Calculates the endpoint of the cable

    pN = get_endpoint(p0, θ, L)
Returns the position of the endpoint.

    pN, vN = get_endpoint(p0, v0, θ, ω, L)
Returns the position and velocity of the endpoint.

Arguments:
  - `p0`: Position of the first link (2D vector)
  - `v0`: Velocity of the first link (2D vector)
  - `θ`: Orientation of the links (`n`D vector, where `n` is the number of segments)
  - `ω`: Angular velocities of the links (`n`D vector, where `n` is the number of segments)
  - `L`: Length of the segments (`n`D vector or scalar; if `L` is scalar, all segments have equal length)
"""
function get_endpoint(p0::Rn, v0::Rn, θ::Rn, ω::Rn, L::Union{Rn,Real})
    xN = p0[1] + sum(L .* cos.(θ))
    yN = p0[2] + sum(L .* sin.(θ))

    v_xN = v0[1] - sum(L .* sin.(θ) .* ω)
    v_yN = v0[2] + sum(L .* cos.(θ) .* ω)

    return [xN; yN], [v_xN; v_yN]
end

function get_endpoint(p0::Rn, θ::Rn, L::Union{Rn,Real})
    xN = p0[1] + sum(L .* cos.(θ))
    yN = p0[2] + sum(L .* sin.(θ))

    return [xN; yN]
end

"""
Calculates the midpoints of the cable

    x, y = get_midpoints(p0, θ, L)
Returns the position of the midpoints.

    x, y, vx, vy = get_midpoints(p0, v0, θ, ω, L)
Returns the position and velocity of the midpoints.

Arguments:
  - `p0`: Position of the first link (2D vector)
  - `v0`: Velocity of the first link (2D vector)
  - `θ`: Orientation of the links (`n`D vector, where `n` is the number of segments)
  - `ω`: Angular velocities of the links (`n`D vector, where `n` is the number of segments)
  - `L`: Length of the segments (`n`D vector or scalar; if `L` is scalar, all segments have equal length)
"""
function get_midpoints(p0::Rn, v0::Rn, θ::Rn, ω::Rn, L::Union{Rn,Real})
    N = length(θ)

    # Let x_kn, y_kn be the knot coordinates.
    # Midpoints are given by:
    #  x_mid = 0.5 * (x_kn[1:end-1] + x_kn[2:end])
    #  y_mid = 0.5 * (y_kn[1:end-1] + y_kn[2:end])
    x, y, v_x, v_y = config2knots(p0, v0, θ, ω, L)
    x .*= 0.5
    x[1:N] .+= x[2:end]
    y .*= 0.5
    y[1:N] .+= y[2:end]
    v_x .*= 0.5
    v_x[1:N] .+= v_x[2:end]
    v_y .*= 0.5
    v_y[1:N] .+= v_y[2:end]

    return x[1:N], y[1:N], v_x[1:N], v_y[1:N]
end

function get_midpoints(p0::Rn, θ::Rn, L::Union{Rn,Real})
    N = length(θ)

    # Let x_kn, y_kn be the knot coordinates.
    # Midpoints are given by:
    #  x_mid = 0.5 * (x_kn[1:end-1] + x_kn[2:end])
    #  y_mid = 0.5 * (y_kn[1:end-1] + y_kn[2:end])
    x, y = config2knots(p0, θ, L)
    x .*= 0.5
    x[1:N] .+= x[2:end]
    y .*= 0.5
    y[1:N] .+= y[2:end]

    return x[1:N], y[1:N]
end
get_midpoints(p0::Rn, θ::Rn, L::Real) = get_midpoints(p0, θ, L .* ones(length(θ)))

"""
Calculates the knots (links) of the cable

    x, y = config2knots(p0, θ, L)
Returns the coordinates of the cable links

    x, y, v_x, v_y = config2knots(p0, v0, θ, ω, L)
Returns the coordinates of the cable links (`x`, `y`) and their velocities (`v_x`, `v_y`).

The coordinates of the links are given by
```
    x[i+1] = x[i] + L[i]*cos(θ[i])
    y[i+1] = y[i] + L[i]*sin(θ[i])

    [x[1]; y[1]] = p0
```

Arguments:
  - `p0`: Position of the first link (2D vector)
  - `v0`: Velocity of the first link (2D vector)
  - `θ`: Orientation of the links (`n`D vector, where `n` is the number of segments)
  - `ω`: Angular velocities of the links (`n`D vector, where `n` is the number of segments)
  - `L`: Length of the segments (`n`D vector or scalar; if `L` is scalar, all segments have equal length)
"""
function config2knots(p0::Rn, v0::Rn, θ::Rn, ω::Rn, L::Union{Rn,Real})
    x = p0[1] .+ [0; cumsum(L.*cos.(θ))]
    y = p0[2] .+ [0; cumsum(L.*sin.(θ))]

    v_x = v0[1] .+ [0; cumsum(-L.*ω.*sin.(θ))]
    v_y = v0[2] .+ [0; cumsum( L.*ω.*cos.(θ))]

    return x, y, v_x, v_y
end

function config2knots(p0::Rn, θ::Rn, L::Union{Rn,Real})
    x = p0[1] .+ [0; cumsum(L.*cos.(θ))]
    y = p0[2] .+ [0; cumsum(L.*sin.(θ))]

    return x, y
end

"""
Converts cable knots to configuration space.

    p0, v0, θ, ω = knots2config(x, y, v_x, v_y)
Returns the position and velocity of the first segment, and the orientation and angular velocities of the cable links.
"""
function knots2config(x::Rn, y::Rn, v_x::Rn, v_y::Rn)
    N = length(x)
    if (length(y) != N) || (length(v_x) != N) || (length(v_y) != N)
        throw(DimensionMismatch("x, y, v_x, and v_y must have the same length"))
    end

    p0 = [x[1]; y[1]]
    v0 = [v_x[1]; v_y[1]]

    dx = diff(x)
    dy = diff(y)
    dvx = diff(v_x)
    dvy = diff(v_y)

    θ = atan.(dy, dx)
    ω = (dx.*dvy - dy.*dvx) ./ (dx.^2 + dy.^2)

    return p0, v0, θ, ω
end
