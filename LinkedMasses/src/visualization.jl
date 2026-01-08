export plot_path!, plot_linked_masses!, mute_color
export plot_simulation, snake_to_capital, make_gif

"""
    plot_path!(ax, path_fcn::Function, s_min::Real, s_max::Real; n_points::Int=100, kwargs...)

Plots a path function.

**Arguments**
- `ax`: Target axes/figure.
- `path_fcn::Function`: A parametric function that defines the path.
- `s_min::Real`: The minimum value of the parameter `s` for the path.
- `s_max::Real`: The maximum value of the parameter `s` for the path.

**Keyword Arguments**
- `n_points::Int=100`: The number of points to sample along the path. Defaults to 100.
- `kwargs...`: Additional keyword arguments passed to the plotting function.
"""
function plot_path!(ax, path_fcn::Function, s_min::Real, s_max::Real; n_points::Int=100, kwargs...)
    # Generate points along the path
    s_values = range(s_min, s_max, length=n_points)
    path_points = hcat((path_fcn(s) for s in s_values)...)

    # Plot the path
    plot!(ax, path_points[2,:], path_points[1,:]; kwargs...)
end

"""
    plot_linked_masses!(ax, state::SimulationState, linked_masses::LinkedMassesParameters; kwargs...)

Plots the linked masses' state on the given axis.

**Arguments**
- `ax`: Target axes/figure.
- `state::SimulationState`: The current state of the simulation.
- `linked_masses::LinkedMassesParameters`: The parameters of the linked masses.
- `kwargs...`: Additional keyword arguments for customizing the plot.

    plot_linked_masses!(ax, x::Rn, linked_masses::LinkedMassesParameters; kwargs...)

Plots the linked masses' state represented by the state vector `x` on the given axis.
"""
function plot_linked_masses!(ax, state::Union{SimulationState,AdaptiveSimulationState}, linked_masses::LinkedMassesParameters; kwargs...)
    # Unpack the state
    θ = state.θ
    p0 = state.p0
    #v0 = state.v0

    # Compute the linked masses' positions
    L = linked_masses.L
    p1 = p0 + L * [cos(θ), sin(θ)]

    # Plot the linked masses
    plot!(ax, [p0[2], p1[2]], [p0[1], p1[1]]; kwargs...)    
end

function plot_linked_masses!(ax, x::Rn, linked_masses::LinkedMassesParameters; kwargs...)
    # Unpack the state
    θ = x[3]
    p0 = x[1:2]

    # Compute the linked masses' positions
    L = linked_masses.L
    p1 = p0 + L * [cos(θ), sin(θ)]

    # Plot the linked masses
    plot!(ax, [p0[2], p1[2]], [p0[1], p1[1]]; kwargs...)    
end

"""
    mute_color(color, amount::Real, min_sat=0.05)

Reduces the saturation of the given `color` by the specified `amount`.
The function converts the input color to HSV, decreases its saturation by `amount` (clamped between 0.05 and 1.0),
and returns the new color in HSV format.

**Arguments**
- `color`: A color object that can be converted to HSV.
- `amount::Real`: The amount by which to decrease the saturation.
- `min_sat::Real=0.05`: The minimum saturation value to clamp to.
"""
function mute_color(color, amount::Real, min_sat::Real=0.05)
    old_hsv = HSV(color)
    s = clamp(old_hsv.s - amount, min_sat, 1.0)
    return HSV(old_hsv.h, s, old_hsv.v)
end

"""
    plot_at_time!(plt, t::Rn, x::Rn, t_i::Real, color, label)

Plots the data `(t, x)` on the given plot `plt`, highlighting the portion up to time `t_i`.

**Arguments**
- `plt`: The plot object to modify.
- `t::Rn`: Array of time values.
- `x::Rn`: Array of data values corresponding to `t`.
- `t_i::Real`: The time up to which the data should be highlighted.
- `color`: The color to use for highlighting.
- `label`: The label for the highlighted portion.

**Description**
Plots the entire `(t, x)` trajectory with a muted version of `color`. Then, highlights the segment from the start up to the first index where `t_i <= t` using the specified `color` and `label`.

**Notes**
- If `t_i` is not found within `t`, only the muted trajectory is plotted.
"""
function plot_at_time!(plt, t::Rn, x::Rn, t_i::Real, color, label)
    plot!(plt, t, x, color=mute_color(color, 0.7), label="")
    i = findfirst(t_i .<= t)
    if !isnothing(i)
        plot!(plt, t[1:i], x[1:i], color=color, label=label)
    end
end

function _path_center(path_fcn::Function, s_min::Real, s_max::Real)
    # Compute the center of the path by averaging the start and end points
    p_start = path_fcn(s_min)
    p_end = path_fcn(s_max)
    return (p_start + p_end) / 2    
end

"""
    plot_simulation(sol, t::Real, parameters::Controller, name="Solution"; kwargs...)

Visualizes the simulation state at a specific time `t`, including the spatial configuration of the linked masses, path-following errors, and angular velocities.

**Arguments**
- `sol`: The solution object returned by the differential equation solver (e.g., `DifferentialEquations.jl`).
- `t::Real`: The specific time instant at which to visualize the system.
- `parameters::Controller`: The controller parameters struct, containing model definitions, path functions, and ocean current information.
- `name::String="Solution"`: The base title for the plot.

**Keyword Arguments**
- `scheme=ColorSchemes.Dark2_8`: The color scheme used for plotting elements.
- `show_time::Bool=true`: If true, appends the current time `t` to the plot title.
- `show_current::Bool=false`: If true, visualizes the ocean current vector.
- `show_path::Bool=true`: If true, draws the reference path on the spatial plot.
- `n_path::Integer=100`: The number of points used to discretize and render the path.
- `vehicle_zoom::Bool=false`: If true, centers the plot on the vehicle (`p0`) and zooms in based on the link length.
- `zoom_range::Real=3.0`: The multiplier for the link length `L` to determine the zoom window size when `vehicle_zoom` is active.

**Returns**
- `plt`: A composite plot object (layout `[a; [b c]]`) containing:
    1.  **Top:** Spatial plot (x vs y) showing the linked masses, reference path, and optional ocean current.
    2.  **Bottom Left:** Path-following errors (e_x, e_y) over time, with a marker at time `t`.
    3.  **Bottom Right:** Link angular velocity over time, with a marker at time `t`.
"""
function plot_simulation(sol, t::Real, parameters::Controller, name="Solution";
    scheme=ColorSchemes.Dark2_8, show_time::Bool=true, show_current::Bool=false, show_path::Bool=true,
    n_path::Integer=100, vehicle_zoom::Bool=false, zoom_range::Real=3.0)
    
    state = unpack_state(sol(t))
    path_fcn = _get_path(parameters)
    linked_masses = _get_model(parameters)
    p_path = path_fcn(state.s)
    s = [u[7] for u in sol.u]

    if show_time 
        title_string = @sprintf "%s (t = %.2f s)" name t
    else
        title_string = name
    end
    plt1 = plot(; title=title_string, aspect_ratio=:equal, xlabel=L"$x$ [m]", ylabel=L"$y$ [m]")
    if show_path
        plot_path!(plt1, path_fcn, minimum(s), maximum(s); color=:black, label="Path", n_points=n_path)
    end
    plot_linked_masses!(plt1, state, linked_masses; marker=:o, color=scheme[1], label="Linked Masses")
    scatter!(plt1, [p_path[2]], [p_path[1]], color=scheme[2], label="Reference")
    if show_current
        p_o = vehicle_zoom ? state.p0 : _path_center(path_fcn, minimum(s), maximum(s))
        V_c = _get_ocean_current(parameters, state.p0, t)
        plot!(plt1, [p_o[2], p_o[2]+V_c[2]], [p_o[1], p_o[1]+V_c[1]], 
            color=scheme[3], label="Ocean Current", arrow=true)
    end
    if vehicle_zoom
        xlims!(plt1, state.p0[2] .+ (-zoom_range*linked_masses.L, zoom_range*linked_masses.L))
        ylims!(plt1, state.p0[1] .+ (-zoom_range*linked_masses.L, zoom_range*linked_masses.L))
    end

    errors = get_errors(sol, parameters)
    plt2 = plot(; title="Path-following error", xlabel=L"$t$ [s]", ylabel=L"$e$ [m]")
    plot_at_time!(plt2, sol.t, errors[1], t, scheme[1], L"e_x")
    plot_at_time!(plt2, sol.t, errors[2], t, scheme[2], L"e_y")

    θ_dot = [u[6] for u in sol.u]
    plt3 = plot(; title="Link angular velocity", xlabel=L"$t$ [s]", ylabel=L"$\dot{\theta}$ [rad/s]")
    plot_at_time!(plt3, sol.t, θ_dot, t, scheme[1], "")

    plt = plot(plt1, plt2, plt3, layout=@layout([a; [b c]]), size=(800, 600))    
    return plt
end

snake_to_capital(name::String) = join(uppercasefirst.(split(name, '_')), " ")

"""
    make_gif(sol, parameters::Controller, name::String; 
             output_dir::String=@__DIR__, fps::Int=20, speed::Real=1.0, kwargs...)

Generates an animated GIF from a simulation solution.

This function iterates through the time steps of the solution `sol`, creating a plot for each frame using `plot_simulation`, and compiles them into a GIF file.

**Arguments**
- `sol`: The solution object returned by the differential equation solver (e.g., `ODESolution`).
- `parameters::Controller`: The controller parameters used during the simulation, passed to the plotting function.
- `name::String`: The base name for the output GIF file (without extension) and the title used in the plot.

**Keyword Arguments**
- `output_dir::String`: The directory where the GIF will be saved. Defaults to the current directory (`@__DIR__`).
- `fps::Int`: Frames per second for the output GIF. Defaults to `20`.
- `speed::Real`: The playback speed multiplier relative to simulation time. Defaults to `1.0` (real-time).
- `kwargs...`: Additional keyword arguments passed directly to `plot_simulation`.

**Returns**
- A `Plots.AnimatedGif` object representing the generated animation.
"""
function make_gif(sol, parameters::Controller, name::String; 
    output_dir::String=@__DIR__, fps::Int=20, speed::Real=1.0, kwargs...)
    

    dt = speed / fps
    t_values = sol.t[1]:dt:sol.t[end]

    anim = @animate for t in t_values
        plot_simulation(sol, t, parameters, snake_to_capital(name); kwargs...)
    end
    return gif(anim, joinpath(output_dir, "$name.gif"), fps=fps)
end
