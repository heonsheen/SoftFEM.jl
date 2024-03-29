using SparseArrays
using LinearAlgebra

# abstract description of an elastic object represented by a mesh
abstract type ElasticObject
### Common Attributes for abstract ElasticObject 
#=
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
    dF::Matrix{Float64} # array of deformation gradient differentials ((dim * NT) x dim)
    W::Vector{Float64} # reference volume of each element (NT x 1)

    T::Matrix{Float64} # mapping vectorized nodal position in a tri to 
                       #   its vectorized deformation gradient (4NT by 6)
                       #   definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)

    M::SparseMatrixCSC{Float64,Int64} # Mass matrix
    K_prev::SparseMatrixCSC{Float64,Int64} # Stiffness matrix of previous timestep
    f_prev::SparseVector{Float64,Int64} # force vector of previous timestep
    K0::SparseMatrixCSC{Float64,Int64} # 

    mat::Material # elastic material description
=#       
end

### methods
function compute_elastic_stiffness_matrix(obj::ElasticObject)
    # abstract definition
end 

function compute_elastic_force(obj::ElasticObject)
    # abstract definition
end

function compute_interface_stiffness_matrix(obj::ElasticObject)
    # abstract definition
end 

function compute_interface_force(obj::ElasticObject)
    # abstract definition
end

function compute_force_differential(obj::ElasticObject)
end

function update_pos(obj::ElasticObject, dx::Vector{Float64})
end