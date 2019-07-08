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

    # damping
    alpha::Float64 # Rayleigh damping multiplier for M
    beta::Float64 # Rayleigh damping multiplier for K

    rho::Float64 # density

    eta::Float64 # DG jump parameter
    DG_IP::Bool # DG IP method flag

### Constructor
    # mp is a Dict containing either YM and PR or Lamer parameters
    function NeohookeanMaterial(
        mp::Dict{String,Float64}, damping::Vector{Float64}, rho::Float64, 
        eta::Float64 = 0.0, DG_IP::Bool = false)
        mp = compute_parameters(mp)
        @assert (size(damping,1) == 2) "damping parameter array must be size 2"

        new(mp["mu"], mp["lambda"], mp["nu"], mp["E"], mp["kappa"], 
            damping[1], damping[2], rho, eta, DG_IP)
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
    dim = size(F,1)

    J = det(F)
    F_inv = inv(F)

    mat.mu * Array{Float64}(I,dim^2,dim^2) + 
        (mat.mu - mat.lambda * log(J)) * kron(F_inv,F_inv) + 
        mat.lambda * vec(F_inv') * vec(F_inv')'
end