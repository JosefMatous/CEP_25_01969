export ASVDimensions, plot_towed_system!
export default_asv_dimensions

"""
    Dimensions of the ASV for visualization purposes.

The ASV is represented as a polygon with 5 vertices: bow, front right, aft right, aft left, front left.

**Fields:**
$(FIELDS)
"""
struct ASVDimensions
    "Body length"
    L_body::Real
    "Bow length"
    L_bow::Real
    "Body width"
    W_body::Real
end

default_asv_dimensions() = ASVDimensions(20.0, 6.0, 8.0)

function asv_vertices(q::Rn, dim::ASVDimensions)
    L_body = dim.L_body
    L_bow = dim.L_bow
    W_body = dim.W_body

    # Vertices of the polygon are in this order:
    #  [bow, front right, aft right, aft left, front left]
    x_center = (3*L_body + L_bow) / 5 # Geometric center of the polygon
    x0 = [L_body + L_bow, L_body, 0, 0, L_body] .- x_center
    y0 = [0, W_body, W_body, -W_body, -W_body] ./ 2

    cpsi = cos(q[3])
    spsi = sin(q[3])

    x = q[1] .+ x0.*cpsi .- y0.*spsi
    y = q[2] .+ x0.*spsi .+ y0.*cpsi

    return x, y
end

const PAYLOAD_MARKER = :pentagon
const CONTROL_POINT_MARKER = :xcross
const PATH_MARKER = :square

function plot_towed_system!(plt::Plots.Plot, pls::PlottingState, dim::ASVDimensions, color)
    x_asv, y_asv = asv_vertices(pls.q_asv, dim)
    Plots.plot!(plt, y_asv, x_asv; seriestype=:shape, color=color, label="", fillopacity=0.5)
    Plots.plot!(plt, pls.y_cable, pls.x_cable, color=color, label="")
    Plots.scatter!(plt, [pls.q_payload[2]], [pls.q_payload[1]], color=color, marker=PAYLOAD_MARKER, label="")
    Plots.scatter!(plt, [pls.p_ε[2]], [pls.p_ε[1]], color=color, marker=CONTROL_POINT_MARKER, label="")
    Plots.scatter!(plt, [pls.p_path[2]], [pls.p_path[1]], color=color, marker=PATH_MARKER, label="")

    return plt
end
