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