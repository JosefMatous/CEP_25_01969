export ASVPoint

const ASVType = Union{ASVModel,ControlledASV}

"""
Represents a point on an ASV, offset a given distance from the center of mass

Fields:
$FIELDS
"""
struct ASVPoint <: EulerLagrangeX.ObjectPoint
    "The ASV"
    obj::EulerLagrangeX.Object
    "Offset from the CO (the point is calculated as: p = [x + d[1]*cos(ψ); y + d[2]*sin(ψ)])"
    d::Rn

    function ASVPoint(obj::EulerLagrangeX.Object, d::Rn)
        if !isa(obj.sys, ASVType) || (EulerLagrangeX.get_num_coordinates(obj.sys) != 3)
            error("The object must be an ASV, i.e., a marine system with three degrees of freedom")
        end
        if length(d) != 2
            throw(DimensionMismatch("d must be a 2D vector"))
        end

        return new(obj, d)
    end
end

EulerLagrangeX.get_num_dimensions(::ASVPoint) = 2

function EulerLagrangeX.get_point(q::Rn, q_dot::Rn, point::ASVPoint)
    R = rotation_3dof_static(q)

    r = q_dot[3]
    S = @SMatrix [0 -r; r 0]

    p = q[1:2] + R * point.d
    v = q_dot[1:2] + R * (S * point.d)

    return p, v
end

function EulerLagrangeX.get_point_acceleration(q::Rn, q_dot::Rn, point::ASVPoint)
    R = rotation_3dof_static(q)
    r = q_dot[3]
    Sd = @SVector [-point.d[2]; point.d[1]]

    # Acceleration:
    #  a = q_ddot[1:2] + R*S^2*d + R*S_dot*d
    # => a0 = R*S^2*d = -R*d.*r^2
    # =>  J = [I(2) R*[-d[2]; d[1]]]

    a0 = (R * point.d) .* (-r^2)
    J = hcat(LinearAlgebra.I(2), R * Sd)

    return J, a0
end
