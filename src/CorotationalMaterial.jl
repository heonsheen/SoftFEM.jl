include("Material.jl")
using LinearAlgebra

# material described by corotational elasticity
mutable struct CorotationalMaterial <: Material
### Attributes   
    mu::Float64 # 2nd Lame parameter (Shear modulus)
    lambda::Float64 # 1st Lame parameter

    # material model
    nu::Float64 # Poisson's ratio
    E::Float64 # Young's modulus
    kappa::Float64 # bulk modulus

    rho::Float64 # density

### Constructor
    # mp is a Dict containing either YM and PR or Lamer parameters
    function CorotationalMaterial(mp::Dict{String,Float64}, rho::Float64)
        mp = compute_parameters(mp)

        new(mp["mu"], mp["lambda"], mp["nu"], mp["E"], mp["kappa"], rho)
    end
end

function compute_energy(F::Matrix{Float64}, mat::CorotationalMaterial)
    # Polar Decomposition
    svd = svd(F)
    R = F.U * F.V'
    S = F.V * F.S * F.V'

    mat.mu * norm(F-R)^2 + mat.lambda / 2 * tr(F.S - I)^2
end

function compute_PK1(F::Matrix{Float64}, mat::CorotationalMaterial)
    # Polar Decomposition
    svd = svd(F)
    R = F.U * F.V'
    S = F.V * F.S * F.V'

    2 * mat.mu * (F-R) + mat.lambda * tr(R'*F - I) * R
end

function compute_dP(F::Matrix{Float64}, dF::Matrix{Float64}, mat::CorotationalMaterial)
    # Polar Decomposition
    svd = svd(F)
    R = F.U * F.V'
    S = F.V * F.S * F.V'

    Eps = [0 1; -1 0] # Levi-Civita symbol
    dR = R * (dot(Eps, inv(tr(S)*I-S)*dot(Eps', R'*dF)))
    2 * mat.mu * dF + mat.lambda * (tr(dR'*F) + tr(R'*dF))*R +
        (mat.lambda * tr(R'*F-I) - 2*mat.mu) * dR
end

function compute_C(F::Matrix{Float64}, mat::CorotationalMaterial)
    # TODO: Need to figure this out...
end