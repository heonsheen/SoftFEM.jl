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
    function LinearElasticMaterial(mp::Dict{String,Float64}, rho::Float64)
        mp = compute_parameters(mp)

        new(mp["mu"], mp["lambda"], mp["nu"], mp["E"], mp["kappa"], rho)
    end
end

function compute_energy(F::Matrix{Float64}, mat::LinearElasticMaterial)
    dim = size(F,1)
    strain = 1/2 * (F + F') - Array{Float64}(I,dim,dim) # infinitesimal strain tensor
    mat.mu * tr(strain'*strain) + mat.lambda/2 * tr(strain)^2
end

function compute_PK1(F::Matrix{Float64}, mat::LinearElasticMaterial)
    dim = size(F,1)
    strain = 1/2 * (F + F') - Array{Float64}(I,dim,dim) #infinitesimal strain tensor
    2 * mat.mu * strain + mat.lambda * tr(strain) * I
end

function compute_dP(F::Matrix{Float64}, dF::Matrix{Float64}, mat::LinearElasticMaterial)
    mat.mu * (dF + dF') + mat.lambda * tr(dF) * I
end

function compute_C(F::Matrix{Float64}, mat::LinearElasticMaterial)
    if size(F,1) == 2
        mat.mu * [2 0 0 0; 0 1 1 0; 0 1 1 0; 0 0 0 2] +
            mat.lambda * [1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 1]
    else
        I3 = Array(I,3,3)
        I9 = Array(I,9,9)
        Kmm =
            [1     0     0     0     0     0     0     0     0;
            0     0     0     1     0     0     0     0     0;
            0     0     0     0     0     0     1     0     0;
            0     1     0     0     0     0     0     0     0;
            0     0     0     0     1     0     0     0     0;
            0     0     0     0     0     0     0     1     0;
            0     0     1     0     0     0     0     0     0;
            0     0     0     0     0     1     0     0     0;
            0     0     0     0     0     0     0     0     1]
        mat.mu * (I9 + Kmm) + mat.lambda * (Kmm * vec(I3) * vec(I3)')
    end
end