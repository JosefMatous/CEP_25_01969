export LinkedMassesParameters, linked_masses_ode
export payload_position

"""
    struct LinkedMassesParameters

Parameters of two masses connected by a rigid link.

**Fields**
- `L::Real`: Length of the link between the two masses.
- `m0::Real`: Mass of the first mass.
- `m::Real`: Mass of the second mass.
- `c0::Real`: Damping coefficient of the first mass.
- `c::Real`: Damping coefficient of the second mass.
"""
struct LinkedMassesParameters
    "Link length"
    L::Real
    "First mass"
    m0::Real
    "Second mass"
    m::Real
    "First damping coefficient"
    c0::Real
    "Second damping coefficient"
    c::Real
end

"""
    linked_masses_ode(x::Rn, f::Rn, params::LinkedMassesParameters, V_c::Rn=zeros(2))

ODEs of two linked masses connected via a rigid link.

**Arguments**
- `x::Rn`: State vector, where
    - `x[1:2]`: Position of the first mass (2D vector).
    - `x[3]`: Angle of the link.
    - `x[4:5]`: Velocity of the first mass (2D vector).
    - `x[6]`: Angular velocity of the link.
- `f::Rn`: External force applied to the first mass (2D vector).
- `params::LinkedMassesParameters`: Structure containing system parameters
- `V_c::Rn=zeros(2)`: Fluid flow velocity (default is zero).

**Returns**
- `dx`: Time derivative of the state vector `x`.
"""
function linked_masses_ode(x::Rn, f::Rn, params::LinkedMassesParameters, V_c::Rn=zeros(2))
    # Unpack parameters
    L = params.L
    m0 = params.m0
    m = params.m
    c0 = params.c0
    c = params.c

    # Unpack states
    θ = x[3]
    v0 = x[4:5]
    θ_dot = x[6]

    # Kinematic equations
    dx = similar(x)
    dx[1:2] = v0 # p0_dot
    dx[3] = θ_dot # θ_dot

    # Dynamic equations
    Γ = [cos(θ), sin(θ)]
    dΓ = [-sin(θ), cos(θ)]
    # Define η = [p0; θ], η_dot = [v0; θ_dot]
    # System equations: M * η_ddot + C * η_dot = J0 * F0 + J * F + J0 * f

    # Mass matrix inverse
    Minv = [(m * dΓ * dΓ' + m0 * I(2)) ./ (m0 * (m + m0))     -dΓ ./ (m0 * L);
            -dΓ' ./ (m0 * L)              (m0 + m) / (m0 * m * L^2)]
    # Coriolis matrix
    C = [m * L * θ_dot^2 * Γ; 0] # = -C * η_dot
    # Hydrodynamic forces
    F0 = -c0 * (v0 - V_c) # Damping force on the first mass
    v = v0 + L * θ_dot * dΓ # Velocity of the second mass
    F = -c * (v - V_c) # Damping force on the second mass
    # Jacobian matrices
    J0 = [I(2); zeros(1, 2)] # ∂p0/∂η
    J = [I(2); L*dΓ'] # ∂p/∂η

    # System equations
    η_ddot = Minv * (C + J0*(F0 + f) + J*F)
    dx[4:6] = η_ddot

    return dx
end

function payload_position(x::Rn, params::LinkedMassesParameters)
    # Unpack parameters
    L = params.L
    p0 = x[1:2]
    θ = x[3]

    # Calculate the position of the payload
    p = p0 + L * [cos(θ), sin(θ)]

    return p    
end
