include("Material.jl")
using LinearAlgebra

# material described by linear elasticity
mutable struct LinearElasticMaterial <: Material
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
    function LinearMaterialModel(mp::Dict{String,Float64}, rho::Float64)
        mp = compute_parameters(mp)

        new(mp["mu"], mp["lambda"], mp["nu"], mp["E"], mp["kappa"], rho)
    end
end

function compute_energy(F::Matrix{Float64}, mat::LinearElasticMaterial)
    strain = 1/2 * (F + F') - Array{Float64}(I,2,2)
    mat.mu * tr(strain'*strain) + mat.lambda/2 * tr(strain)^2
end

function compute_PK1(F::Matrix{Float64}, mat::LinearElasticMaterial)
    strain = 1/2 * (F + F') - Array{Float64}(I,2,2)
    2 * mat.mu * strain + mat.lambda * tr(strain) * I
end