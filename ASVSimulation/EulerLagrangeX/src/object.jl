export Object, ObjectPoint
export get_parent, get_type, is_on, get_num_dimensions, get_point, get_point_acceleration


const ObjType = Union{EulerLagrangeSystem,ControlledEulerLagrangeSystem}

"""
A unique object

A structure consisting of an Euler-Lagrange system and a unique ID number.
"""
struct Object{T}
    sys::T
    id::Integer

    """
    Creates a new object

        obj = Object(sys, id)
    Creates an object out of the given Euler-Lagrange system.
    NB: The function does not check if `id` is unique.
    """
    function Object(sys::T, id::Integer) where T
        if !isa(sys, ObjType)
            throw(TypeError(:Object, ObjType, T))
        end
        return new{T}(sys, id)
    end
end

"Abstract type representing a point on an Euler-Lagrange system"
abstract type ObjectPoint end

"""
Returns the parent of the point.

    obj = get_parent(point)
Returns the parent object that the point belongs to.
NB: This function must be implemented in the derived types.
"""
function get_parent(point::T) where T<:ObjectPoint
    if hasfield(T, :obj)
        return point.obj
    else
        error("get_parent must be implemented for the derived type")
    end
end

"""
Returns the type of the parent object

    T = get_type(point)
Returns the type `T <: Union{EulerLagrangeSystem,ControlledEulerLagrangeSystem}`
of the parent object
"""
get_type(point::ObjectPoint) = typeof(get_parent(point).sys)

"""
Checks if the given point is on the given object

    flag = is_on(point, obj)
Checks if `point::ObjectPoint` is the child of `obj::Object`
"""
is_on(point::ObjectPoint, obj::Object) = (get_parent(point) === obj)

"""
Returns the number of dimensions of the point
"""
function get_num_dimensions(::ObjectPoint)
    error("get_num_dimensions must be implemented for the derived type")
end

"""
Calculates the position and velocity of the point

    p, v = get_point(q, q_dot, point)
Calculates the position (`p`) and velocity (`v`) of the `point`, given 
generalized coordinates `q` and generalized velocities `q_dot`.
"""
function get_point(::Rn, ::Rn, ::ObjectPoint)
    error("get_point must be implemented for the derived type")
end

"""
Calculates the acceleration of the point

    J, a0 = get_point_acceleration(q, q_dot, point)
Calculates the matrix `J` and vector `a0` such that the acceleration of the point
is given by: `a0 + J*q_ddot`, where `q_ddot` is the acceleration of the system.
"""
function get_point_acceleration(q::Rn, q_dot::Rn, point::ObjectPoint)
    vel_q_fcn = (x) -> get_point(x, q_dot, point)[2]
    vel_q_dot_fcn = (x) -> get_point(q, x, point)[2]

    # Let v(q, q_dot) be the velocity of the point.
    # Then, the acceleration of the point is:
    #  a = ∂v/∂q * q_dot  +  ∂v∂q_dot * q_ddot
    # => a0 = ∂v/∂q * q_dot
    # =>  J = ∂v∂q_dot

    # Calculate a0
    a0 = ForwardDiff.jacobian(vel_q_fcn, q) * q_dot
    # Calculate J
    J = ForwardDiff.jacobian(vel_q_dot_fcn, q_dot) # ∂v∂q_dot
    return J, a0
end
