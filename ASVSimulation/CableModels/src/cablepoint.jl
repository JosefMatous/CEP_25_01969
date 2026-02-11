export CablePoint
export cable_begin, cable_end

"""
Represents a point on a rigid cable

Fields:
$(FIELDS)
"""
struct CablePoint <: EulerLagrangeX.ObjectPoint
    "The cable"
    obj::EulerLagrangeX.Object{Cable}
    "Index of the cable segment (an integer between 1 and `N`, where `N` is the number of segments)"
    segment::Integer
    "Relative distance (a real number between 0 and 1; 0 means the beginning of the segment, 1 means the end of the segment)"
    d_rel::Real

    function CablePoint(obj::EulerLagrangeX.Object{Cable}, segment::Integer, d_rel::Real)
        N = obj.sys.N
        if (segment < 1) || (segment > N)
            throw(DomainError(segment, "segment must be an integer between 1 and $(N)"))
        end
        if (d_rel < 0) || (d_rel > 1)
            throw(DomainError(d_rel, "d_rel must be a real number between 0 and 1"))
        end

        return new(obj, segment, d_rel)
    end
end

EulerLagrangeX.get_num_dimensions(::CablePoint) = 2

function EulerLagrangeX.get_point(q::Rn, q_dot::Rn, point::CablePoint)
    sys = point.obj.sys
    L = sys.L
    vector_L = isa(L, Rn)
    #N_q = EulerLagrangeX.get_num_coordinates(sys)

    # Let:
    #  q = [p0; θ]
    #  q_dot = [v0; ω]
    # Then, the position of the point is:
    #  p = p0 + ∑_{j=1}^{segment-1}L_j*[cos(θ[j]), sin(θ[j])] + d_rel*L_s*[cos(θ[s]), sin(θ[s])]
    # and the velocity is:
    #  v = v0 + ∑_{j=1}^{segment-1}L_j*ω_j*[-sin(θ[j]), cos(θ[j])] + d_rel*L_s*ω[s]*[-sin(θ[s]), cos(θ[s])]

    p = q[1:2] # p0
    v = q_dot[1:2] # v0

    # Add the sums
    for j = 1:(point.segment-1)
        L_j = vector_L ? L[j] : L

        cθ = cos(q[j+2]); sθ = sin(q[j+2])
        ω = q_dot[j+2]

        p[1] += L_j*cθ
        p[2] += L_j*sθ

        v[1] -= L_j*ω*sθ
        v[2] += L_j*ω*cθ
    end

    # Add the last segment
    j = point.segment
    L_j = point.d_rel * (vector_L ? L[j] : L)

    cθ = cos(q[j+2]); sθ = sin(q[j+2])
    ω = q_dot[j+2]

    p[1] += L_j*cθ
    p[2] += L_j*sθ

    v[1] -= L_j*ω*sθ
    v[2] += L_j*ω*cθ

    return p, v
end

function EulerLagrangeX.get_point_acceleration(q::Rn, q_dot::Rn, point::CablePoint)
    sys = point.obj.sys
    L = sys.L
    vector_L = isa(L, Rn)
    N_q = length(q)

    # Let:
    #  q = [p0; θ]
    #  q_dot = [v0; ω]
    # Then, acceleration of the point is:
    #  a = dv0/dt + ∑_{j=1}^{segment-1}L_j*ω_j^2*[-cos(θ[j]), -sin(θ[j])] + d_rel*L_s*ω[s]^2*[-cos(θ[s]), -sin(θ[s])]
    #             + ∑_{j=1}^{segment-1}L_j*dω_j/dt*[-sin(θ[j]), cos(θ[j])] + d_rel*L_s*dω[s]/dt*[-sin(θ[s]), cos(θ[s])]
    #
    # => a0 = ∑_{j=1}^{segment-1}L_j*ω_j^2*[-cos(θ[j]), -sin(θ[j])] + d_rel*L_s*ω[s]^2*[-cos(θ[s]), -sin(θ[s])]
    # =>  J = [1 0 -L_1*sin(θ[1]) -L_2*sin(θ[2]) ... -d_rel*L_s*sin(θ[s]) 0 ... 0
    #          0 1  L_1*cos(θ[1])  L_2*cos(θ[2]) ...  d_rel*L_s*cos(θ[s]) 0 ... 0]

    # Initialize
    T = promote_type(eltype(q), eltype(q_dot))
    J = zeros(T, 2, N_q)
    a0 = zeros(T, 2)

    J[1,1] = 1
    J[2,2] = 1

    # Add the sums
    for j = 1:(point.segment-1)
        L_j = vector_L ? L[j] : L

        cθ = cos(q[j+2]); sθ = sin(q[j+2])
        ω2 = q_dot[j+2]^2

        a0[1] -= L_j*ω2*cθ
        a0[2] -= L_j*ω2*sθ

        J[1,j+2] = -L_j*sθ
        J[2,j+2] =  L_j*cθ
    end

    # Add the last segment
    j = point.segment
    L_j = point.d_rel * (vector_L ? L[j] : L)

    cθ = cos(q[j+2]); sθ = sin(q[j+2])
    ω2 = q_dot[j+2]^2

    a0[1] -= L_j*ω2*cθ
    a0[2] -= L_j*ω2*sθ

    J[1,j+2] = -L_j*sθ
    J[2,j+2] =  L_j*cθ

    return J, a0
end

"Creates a CablePoint representing the beginning of the cable"
function cable_begin(obj::EulerLagrangeX.Object{Cable})
    return CablePoint(obj, 1, 0.)
end

"Creates a CablePoint representing the end of the cable"
function cable_end(obj::EulerLagrangeX.Object{Cable})
    N = obj.sys.N
    return CablePoint(obj, N, 1.)
end
