export ASVModel

"""
Parameters of an ASV model with linear and quadratic hydrodynamic damping.

$(FIELDS)
"""
struct ASVModel <: EulerLagrangeX.EulerLagrangeSystem
    "ASV mass"
    m::Real   
    "Yaw inertia"
    J::Real   
    "Position of the center of mass in the body-fixed coordinate frame"
    xG::Real  
    "Added mass in the surge direction"
    Xu::Real  
    "Added mass in the sway direction"
    Yv::Real  
    "Added mass in yaw"
    Nr::Real  
    "Cross-coupling between sway and yaw"
    Yr::Real  
    "Linear hydrodynamic damping in surge"
    Dul::Real 
    "Linear hydrodynamic damping in sway"
    Dvl::Real 
    "Linear hydrodynamic damping in yaw"
    Drl::Real 
    "Quadratic hydrodynamic damping in surge"
    Duq::Real 
    "Quadratic hydrodynamic damping in sway"
    Dvq::Real 
    "Quadratic hydrodynamic damping in yaw"
    Drq::Real 
    "Ocean current"
    V::Union{Rn, Function}
    "Perturbation"
    μ::Union{Rn, Function}
end

EulerLagrangeX.get_num_coordinates(::ASVModel) = 3

"""
Returns the inertia and added mass matrices
"""
function inertia_matrices(mdl::ASVModel)
    # Rigid-body
    M = @SMatrix [
        mdl.m            0                    0;
            0        mdl.m         mdl.m*mdl.xG;
            0 mdl.m*mdl.xG mdl.m*mdl.xG^2+mdl.J
    ]
    # Added
    MA = @SMatrix [
        mdl.Xu      0      0;
             0 mdl.Yv mdl.Yr;
             0 mdl.Yr mdl.Nr
    ]

    return M, MA
end

function EulerLagrangeX.lagrangian(q::Rn, q_dot::Rn, t::Real, mdl::ASVModel)
    M, MA = inertia_matrices(mdl)

    J = transform_3dof_static(q)
    ν = J'*q_dot # body-fixed velocities
    ν_rel = ν - J'*[ocean_current(t, mdl.V); 0] # relative velocities (w.r.t. Ocean current)

    ℓ = 0.5*(ν'*M*ν + ν_rel'*MA*ν_rel)
    return ℓ
end

function EulerLagrangeX.forces(q::Rn, q_dot::Rn, t::Real, mdl::ASVModel)
    # Hydrodynamic forces
    J = transform_3dof_static(q)
    ν_rel = J' * (q_dot - [ocean_current(t, mdl.V);0]) # Relative velocity in the body-fixed frame

    DL = @SVector [mdl.Dul; mdl.Dvl; mdl.Drl] # Linear damping coefficients
    DQ = @SVector [mdl.Duq; mdl.Dvq; mdl.Drq] # Quadratic damping coefficients

    D_b = -ν_rel .* (DL + DQ.*abs.(ν_rel)) # Hydrodynamic damping in the body-fixed frame
    D = J * D_b # Hydrodynamic damping in the inertial frame

    # Perturbation
    p = q[1:2] # Position of the ASV
    return D + [mdl.m * perturbation(t, p, mdl.μ); 0]
end
