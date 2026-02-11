export LinkedMassesParameters

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
