"""3DOF transformation matrix; return a static matrix"""
function transform_3dof_static(η::Rn)
    ψ = η[3]
    cψ = cos(ψ); sψ = sin(ψ)

    J = @SMatrix [
        cψ -sψ 0;
        sψ  cψ 0;
         0   0 1
    ]
    return J
end

"""3DOF rotation matrix; return a static matrix"""
function rotation_3dof_static(η::Rn)
    ψ = η[3]
    cψ = cos(ψ); sψ = sin(ψ)

    R = @SMatrix [
        cψ -sψ;
        sψ  cψ
    ]

    return R
end
