# abstract description of an elastic object represented by a mesh
mutable struct ElasticObject
### Attributes    
    N::Int64 # number of nodes
    NT::Int64 # number of elements
    dim::Int64 # dimension of the object
    
    x_node::Matrix{Float64} # nodal position in world coordinates (N x dim)
    X_node::Matrix{Float64} # nodal position in material coordinates (N x dim)
    ec::Matrix{Int64} # element connectivity (NT x {3,4,6,...})

    x::Vector{Float64} # current position in world coordinates ((dim * N) x 1)
    X::Vector{Float64} # reference position in world coordinates ((dim * N) x 1)
    v::Vector{Float64} # velocity in world coordinates ((dim * N) x 1)

    Ds::Matrix{Float64} # array of deformed shape matrix ((dim * NT) x dim)
    Dm::Matrix{Float64} # array of reference shape matrix ((dim * NT) x dim)
    Dm_inv::Matrix{Float64} # array of Dm inverse ((dim * NT) x dim)
    F::Matrix{Float64} # array of deformation gradient ((dim * NT) x dim)
    F_inv::Matrix{Float64} # array of inverse deformation gradient ((dim * NT) x dim)
    W::Vector{Float64} # reference volume of each element (NT x 1)

    M::Matrix{Float64} # Mass matrix
    K_prev::Matrix{Float64} # Stiffness matrix of previous timestep
    f_prev::Matrix{Float64} # force vector of previous timestep
    K0::Matrix{Float64} # 
end
