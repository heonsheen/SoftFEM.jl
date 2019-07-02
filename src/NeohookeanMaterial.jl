include("Material.jl")
using LinearAlgebra

# material described by Neohookean elasticity
mutable struct NeohookeanMaterial <: Material
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
    function NeohookeanMaterial(mp::Dict{String,Float64}, rho::Float64)
        mp = compute_parameters(mp)

        new(mp["mu"], mp["lambda"], mp["nu"], mp["E"], mp["kappa"], rho)
    end
end

function compute_energy(F::Matrix{Float64}, mat::NeohookeanMaterial)
    J = det(F)

    mat.mu / 2 * (tr(F'*F) - 3) - mat.mu * log(J) + mat.lambda / 2 * log(J)^2
end

function compute_PK1(F::Matrix{Float64}, mat::NeohookeanMaterial)
    J = det(F)
    F_inv = inv(F)

    mat.mu * (F - F_inv') + mat.lambda * log(J) * F_inv'
end

function compute_dP(F::Matrix{Float64}, dF::Matrix{Float64}, mat::NeohookeanMaterial)
    J = det(F)
    F_inv = inv(F)
    FidFFi = F_inv * dF * F_inv

    mat.mu * (dF + FidFFi') + mat.lambda * (tr(FidFFi)*F_inv' - log(J)*FidFFi')
end

function compute_C(F::Matrix{Float64}, mat::NeohookeanMaterial)
    J = det(F)
    F_inv = inv(F)

    mat.mu * Array{Float64}(I,4,4) + (mat.mu - mat.lambda * log(J)) * kron(F_inv,F_inv)
        + mat.lambda * vec(F_inv') * vec(F_inv')'
end