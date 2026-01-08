export SimulationState
export AdaptiveSimulationState
export unpack_state, pack_state

"""
    struct SimulationState

Current state of the simulation.

**Fields**
$(FIELDS)
"""
struct SimulationState
    "Position of the first mass"
    p0::Rn
    "Angle of the link"
    θ::Real
    "Velocity of the first mass"
    v0::Rn
    "Angular velocity of the link"
    θ_dot::Real
    "Path parameter"
    s::Real
end

"""
    AdaptiveSimulationState

Current state of the simulation with a model reference adaptive controller (MRAC).

**Fields**
$(FIELDS)
"""
struct AdaptiveSimulationState
    "Position of the first mass"
    p0::Rn
    "Angle of the link"
    θ::Real
    "Velocity of the first mass"
    v0::Rn
    "Angular velocity of the link"
    θ_dot::Real
    "Path parameter"
    s::Real
    "Parameter estimates"
    ζ_hat::Rn
end


"""
    pack_state(st::SimulationState) -> Vector

Packs the current state of the simulation into a vector.

**Arguments**
- `st::SimulationState`: The current state of the simulation to be packed.

**Returns**
- A vector representation of the simulation state.

    pack_state(st::AdaptiveSimulationState) -> Vector

Packs the current state of the simulation with a model reference adaptive controller (MRAC) into a vector.

**Arguments**
- `st::AdaptiveSimulationState`: The current state of the simulation to be packed.

**Returns**
- A vector representation of the simulation state. Elements are ordered as follows:

**See Also**
`unpack_state`
"""
function pack_state(st::SimulationState)
    # Pack the state into a single vector
    return [st.p0; st.θ; st.v0; st.θ_dot; st.s]    
end

function pack_state(st::AdaptiveSimulationState)
    # Pack the state into a single vector
    return [st.p0; st.θ; st.v0; st.θ_dot; st.s; st.ζ_hat]        
end

"""
    unpack_state(x::Rn)

Unpacks the state vector `x` into a `SimulationState`. The inverse of `pack_state`.

**Arguments**
- `x::Rn`: The state vector.

**Returns**
A `SimulationState` or an `AdaptiveSimulationState` object, depending on the length of `x`.

**See Also**
`pack_state`
"""
function unpack_state(x::Rn)
    p0 = x[1:2]
    θ = x[3]
    v0 = x[4:5]
    θ_dot = x[6]
    s = x[7]
    if length(x) > 7
        ζ_hat = x[8:end]
        return AdaptiveSimulationState(p0, θ, v0, θ_dot, s, ζ_hat)
    else
        return SimulationState(p0, θ, v0, θ_dot, s)
    end
end
