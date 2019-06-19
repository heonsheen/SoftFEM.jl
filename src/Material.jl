# abstract description of elastic Material
abstract type Material
### Common Attributes for abstract Material 
#=
    # Lame parameters
    mu::Float64 # 2nd Lame parameter (Shear modulus)
    lambda::Float64 # 1st Lame parameter

    # material model
    nu::Float64 # Poisson's ratio
    E::Float64 # Young's modulus
    kappa::Float64 # bulk modulus

    rho::Float64 # density
=#
end

### methods

# Compute all material parameters from 
#   1. Young's Modulus and Poisson's Ratio
#   2. Lame parameters
function compute_parameters(mp::Dict{String,Float64})
    if get(mp, "E", -1) != -1 && get(mp, "nu", -1) != -1
        E = mp["E"]
        nu = mp["nu"]
        mp["mu"] = 0.5 * E / (1 + nu)
        mp["lambda"] = E * nu / ((1 + nu) * (1 - 2 * nu))
        mp["kappa"] = E / (3 * (1 - 2 * nu))
    elseif get(mp, "mu", -1) != -1 && get(mp, "lambda", -1) != -1
        mu = mp["mu"]
        lambda = mp["lambda"]
        mp["E"] = mu * (3 * lambda + 2 * mu) / (lambda + mu)
        mp["nu"] = lambda / (2 * (lambda + mu))
        mp["kappa"] = lambda + 2 * E / 3
    end
    
    mp
end
