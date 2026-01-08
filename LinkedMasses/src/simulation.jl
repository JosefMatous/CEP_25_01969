export simulate, get_errors


"""
    simulate(ctrl::Controller, x0::Rn, T_stop::Real; solver=Tsit5(), reltol=1e-3)
    simulate(ctrl::Controller, x0::Union{SimulationState,AdaptiveSimulationState}, T_stop::Real; kwargs...)

Simulate the closed-loop system dynamics for a given controller and initial state.

**Arguments**
- `ctrl::Controller`: The controller object defining the control law.
- `x0`: The initial state of the system. This can be provided as:
    - A raw numerical vector (`Rn`).
    - A structured state object (`SimulationState` or `AdaptiveSimulationState`), which will be packed into a vector automatically.
- `T_stop::Real`: The final simulation time (simulation runs from `t=0.0` to `t=T_stop`).

# Keyword Arguments
- `solver`: The differential equation solver to use (default: `Tsit5()`).
- `reltol`: The relative tolerance for the solver (default: `1e-3`).

**Returns**
- `sol`: The solution object.
"""
function simulate(ctrl::Controller, x0::Rn, T_stop::Real; solver=Tsit5(), kwargs...)
    prob = ODEProblem(closed_loop_ode!, x0, (0.0, T_stop), ctrl)
    sol = solve(prob, solver; kwargs...)
    return sol
end

function simulate(ctrl::Controller, x0::Union{SimulationState,AdaptiveSimulationState}, T_stop::Real; kwargs...)
    return simulate(ctrl, pack_state(x0), T_stop; kwargs...)
end

"""
    get_errors(sol, params)

Returns the path-following and velocity errors of the ODE solution.

**Arguments**
- `sol::ODESolution`: ODE solution.
- `params::Controller`: Controller parameters.

**Returns**
- `e_x::Vector{Float64}`: Path-following error in the x-direction.
- `e_y::Vector{Float64}`: Path-following error in the y-direction.
- `v_x::Vector{Float64}`: Velocity error in the x-direction.
- `v_y::Vector{Float64}`: Velocity error in the y-direction.
"""
function get_errors(sol, params::Controller)
    e_x = Float64[]
    e_y = Float64[]
    v_x = Float64[]
    v_y = Float64[]
    path_fcn = _get_path(params)
    los = _get_los(params)
    for u in sol.u
        simstate = unpack_state(u)
        p, v = control_point(simstate, params)
        p_path = path_fcn(simstate.s)
        ∂p_path = ForwardDiff.derivative(path_fcn, simstate.s)
        ψ = atanv(∂p_path)
        R = [cos(ψ) -sin(ψ); sin(ψ) cos(ψ)]
        e = R' * (p - p_path)
        push!(e_x, e[1])
        push!(e_y, e[2])

        v_LOS, _ = line_of_sight(p, simstate.s, path_fcn, los)
        push!(v_x, v[1] - v_LOS[1])
        push!(v_y, v[2] - v_LOS[2])
    end

    return e_x, e_y, v_x, v_y
end

function control_point(simstate::Union{SimulationState,AdaptiveSimulationState}, ctrl::Union{LinearizingController,AdaptiveLinearizingController})
    θ = simstate.θ
    p0 = simstate.p0
    v0 = simstate.v0
    θ_dot = simstate.θ_dot

    lm = ctrl.model
    ε = ctrl.ε

    Γ = [cos(θ), sin(θ)]
    dΓ = [-sin(θ), cos(θ)]
    p = p0 + ε * lm.L * Γ
    v = v0 + ε * lm.L * θ_dot * dΓ
    return p, v    
end
