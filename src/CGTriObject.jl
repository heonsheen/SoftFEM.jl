include("ElasticObject.jl")
include("Geometry.jl")
include("Material.jl")
include("LinearElasticMaterial.jl")

mutable struct CGTriObject <: ElasticObject
### Attributes 
# TODO: How do I make this inherit all attributes from ElasticObject?
# (maybe Lazy::@forward?)
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

    M::SparseMatrixCSC{Float64,Int64} # Mass matrix
    K_prev::SparseMatrixCSC{Float64,Int64} # Stiffness matrix of previous timestep
    f_prev::SparseMatrixCSC{Float64,Int64} # force vector of previous timestep
    K0::SparseMatrixCSC{Float64,Int64} # 

    mat::Material # elastic material description
    
### Constructor
    function CGTriObject(
        mesh::Mesh,
        material::Material)
        N = mesh.n_vertices
        NT = mesh.n_faces
        dim = mesh.dim

        x_node = [vtx.x[i] for vtx in mesh.vertices, i in 1:dim]
        X_node = x_node
        ec = mesh.ec

        x = reshape(x_node, (dim*N, 1))
        X = x_node
        v = zeros((dim*N, 1))

        G = [1 0; 0 1; -1 -1]

        M_I = M_J = Vector{Int64}()
        M_V = Vector{Float64}()
        Kp_I = Kp_J = Vector{Int64}()
        Kp_V = Vector{Float64}()
        fp_I = Vector{Int64}()
        fp_V = Vector{Float64}()
        K0_I = K0_J = Vector{Int64}()
        K0_V = Vector{Float64}()

        for t in 1:NT
            X_t = X_node[[ec[t,i] for i in 1:3],:]
            Dm[2*t-1:2*t,:] = X_t' * G
            Ds[2*t-1:2*t,:] = X_t' * G
            Dm_inv[2*t-1:2*t,:] = inv(X_t' * G)
            F[2*t-1:2*t,:] = Matrix{Float64}(I,2,2)
            F_inv[2*t-1:2*t,:] = Matrix{Float64}(I,2,2)
            vol = det(X_t' * G)/2; # undeformed volume from matrix determinant
            @assert vol > 0 "Wrong orientation in mesh"
            W[t] = vol;

            # Mass Matrix 
            # 1. regular FEM mass matrix
            #=
            Mh_T = abs(vol) / 12 * [2 1 1; 1 2 1; 1 1 2] * rho
            for i = 1:3, j = 1:3
                M_I.push(ec[t,i])
                M_J.push(ec[t,j])
                M_V.push(Mh_T[i,j])
            end
            =# 
            # 2. Lumped mass matrix
            Mh_T = abs(vol) / 12 * [4 0 0; 0 4 0; 0 0 4] * material.rho
            for i = 1:3, j = 1:3
                M_I.push(ec[t,i])
                M_J.push(ec[t,j])
                M_V.push(Mh_T[i,j])
            end


        end
    end
end

function compute_stiffness_matrix(obj::CGTriObject)

end 

function compute_elastic_force(obj::CGTriObject)
    f = zeros(2*obj.N,1)
    for t in 1:obj.NT
        F_t = obj.F[2*t-1:2*t,:]
        Dm_inv_t = obj.Dm_inv[2*t-1:2*t,:]
        H = -obj.W[t] * compute_PK1(F_t, obj.mat) * Dm_inv_t'
        f1 = H[:,1]
        f2 = H[:,2]
        f3 = -(f1 + f2)
        T = obj.ec[t,:]
        f[2*T[1]-1:2*T[1]] += f1
        f[2*T[2]-1:2*T[2]] += f2
        f[2*T[3]-1:2*T[3]] += f3
    end
end

function compute_force_differential(obj::CGTriObject)
    df = zeros(2*obj.N,1)
    for t in 1:obj.NT
        dF_t = obj.dF[2*t-1:2*t,:]
        Dm_inv_t = obj.Dm_inv[2*t-1:2*t,:]
        H = -obj.W[t] * compute_dP(dF_t, obj.mat) * Dm_inv_t'
        df1 = H[:,1]
        df2 = H[:,2]
        df3 = -(df1 + df2)
        T = obj.ec[t,:]
        df[2*T[1]-1:2*T[1]] += df1
        df[2*T[2]-1:2*T[2]] += df2
        df[2*T[3]-1:2*T[3]] += df3
    end
end    
