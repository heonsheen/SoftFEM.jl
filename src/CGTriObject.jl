include("ElasticObject.jl")
include("Geometry.jl")
include("Material.jl")

mutable struct CGTriObject <: ElasticObject
### Attributes 
# TODO: How do I make this automatically inherit all attributes from ElasticObject?
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
    dF::Matrix{Float64} # array of deformation gradient differentials ((dim * NT) x dim)
    W::Vector{Float64} # reference volume of each element (NT x 1)

    T::Matrix{Float64} # mapping vectorized nodal position in a tri to 
                       #   its vectorized deformation gradient (4NT by 6)
                       #   definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)

    M::SparseMatrixCSC{Float64,Int64} # Mass matrix
    K_prev::SparseMatrixCSC{Float64,Int64} # Stiffness matrix of previous timestep
    f_prev::SparseMatrixCSC{Float64,Int64} # force vector of previous timestep
    K0::SparseMatrixCSC{Float64,Int64} # 

    mat::Material # elastic material description
    
### Constructor
    function CGTriObject(
        mesh::Mesh,
        mat::Material)
        N = mesh.n_vertices
        NT = mesh.n_faces
        dim = mesh.dim

        x_node = [vtx.x[i] for vtx in mesh.vertices, i in 1:dim]
        X_node = x_node
        ec = mesh.ec

        x = vec(reshape(x_node', (dim*N, 1)))
        X = vec(reshape(x_node', (dim*N, 1)))
        v = zeros(dim*N)

        G = [1 0; 0 1; -1 -1]
        I2 = Matrix{Float64}(I,2,2)

        M_I = Vector{Int64}()
        M_J = Vector{Int64}()
        M_V = Vector{Float64}()
        K0_I = zeros(Int64, 9*4*NT)
        K0_J = zeros(Int64, 9*4*NT)
        K0_V = zeros(Float64, 9*4*NT)

        Ds = zeros(dim*NT, dim)
        Dm = zeros(dim*NT, dim)
        Dm_inv = zeros(dim*NT, dim)
        F = zeros(dim*NT, dim)
        F_inv = zeros(dim*NT, dim)
        dF = zeros(dim*NT, dim)
        W = zeros(NT)
        
        T = zeros(4 * NT, 6)

        nnz = 0 

        for t in 1:NT
            X_t = X_node[[ec[t,i] for i in 1:3],:]
            Dm[2*t-1:2*t,:] = X_t' * G
            Ds[2*t-1:2*t,:] = X_t' * G
            Dm_inv[2*t-1:2*t,:] = inv(X_t' * G)
            F[2*t-1:2*t,:] = I2
            F_inv[2*t-1:2*t,:] = I2
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
            Mh_T = abs(vol) / 12 * [4 0 0; 0 4 0; 0 0 4] * mat.rho
            for i = 1:3, k = 1:2
                push!(M_I, 2*(ec[t,i]-1) + k)
                push!(M_J, 2*(ec[t,i]-1) + k)
                push!(M_V, Mh_T[i,i])
            end

            T[4*(t-1)+1:4*t,:] = kron((G * inv(X_t' * G))', Matrix{Float64}(I,2,2))

            T_t = T[4*(t-1)+1:4*t,:]
            C = compute_C(I2, mat)

            # Stiffness Matrix
            K_t = W[t] * T_t' * C * T_t
            K_t = 1/2 * (K_t + K_t')

            for i in 1:3, j in 1:3
                K0_I[nnz+1:nnz+4] = repeat((2*(ec[t,i]-1)+1:2*ec[t,i]), 2, 1)
                K0_J[nnz+1:nnz+4] = reshape(repeat((2*(ec[t,j]-1)+1:2*ec[t,j])', 2, 1), 4, 1)
                K0_V[nnz+1:nnz+4] = reshape(K_t[2*(i-1)+1:2*i,2*(j-1)+1:2*j],4,1)
                nnz += 4
            end
        end

        M = sparse(M_I, M_J, M_V, 2*N, 2*N)
        K0 = sparse(K0_I, K0_J, K0_V, 2*N, 2*N)
        K_prev = spzeros(2*N, 2*N)
        f_prev = spzeros(2*N, 1)
        
        new(N, NT, dim, x_node, X_node, ec, x, X, v, 
            Ds, Dm, Dm_inv, F, F_inv, dF, W, T, M, K_prev, f_prev, K0, mat)
    end
end

function compute_stiffness_matrix(obj::CGTriObject)
    I = zeros(Int64,9*4*obj.NT)
    J = zeros(Int64,9*4*obj.NT)
    V = zeros(Float64,9*4*obj.NT)
    nnz = 0
    for t in 1:obj.NT
        F_t = obj.F[2*t-1:2*t,:]
        T_t = obj.T[4*(t-1)+1:4*t,:]
        C = compute_C(F_t, obj.mat)

        K_t = obj.W[t] * T_t' * C * T_t
        K_t = 1/2 * (K_t + K_t')

        for i in 1:3, j in 1:3
            I[nnz+1:nnz+4] = repeat((2*(obj.ec[t,i]-1)+1:2*obj.ec[t,i]), 2, 1)
            J[nnz+1:nnz+4] = reshape(repeat((2*(obj.ec[t,j]-1)+1:2*obj.ec[t,j])', 2, 1), 4, 1)
            V[nnz+1:nnz+4] = reshape(K_t[2*(i-1)+1:2*i,2*(j-1)+1:2*j],4,1)
            nnz += 4
        end
    end

    sparse(I, J, V, 2*obj.N, 2*obj.N)
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

    f
end

function compute_force_differential(obj::CGTriObject)
    df = zeros(2*obj.N,1)
    for t in 1:obj.NT
        F_t = obj.F[2*t-1:2*t,:]
        dF_t = obj.dF[2*t-1:2*t,:]
        Dm_inv_t = obj.Dm_inv[2*t-1:2*t,:]
        H = -obj.W[t] * compute_dP(F_t, dF_t, obj.mat) * Dm_inv_t'
        df1 = H[:,1]
        df2 = H[:,2]
        df3 = -(df1 + df2)
        T = obj.ec[t,:]
        df[2*T[1]-1:2*T[1]] += df1
        df[2*T[2]-1:2*T[2]] += df2
        df[2*T[3]-1:2*T[3]] += df3
    end

    df
end    

function update_pos(obj::CGTriObject, dx::Vector{Float64})
    obj.x = obj.X + dx
    obj.x_node = reshape(obj.x, (obj.dim, obj.N))'
    dx_node = reshape(dx, (obj.dim, obj.N))'

    G = [1 0; 0 1; -1 -1]

    for t in 1:obj.NT
        x_t = obj.x_node[[obj.ec[t,i] for i in 1:3],:]
        dx_t = dx_node[[obj.ec[t,i] for i in 1:3],:]
        obj.Ds[2*t-1:2*t,:] = x_t' * G
        obj.F[2*t-1:2*t,:] = obj.Ds[2*t-1:2*t,:] * obj.Dm_inv[2*t-1:2*t,:]
        obj.F_inv[2*t-1:2*t,:] = obj.Dm[2*t-1:2*t,:] * inv(obj.Ds[2*t-1:2*t,:])
        obj.dF[2*t-1:2*t,:] = dx_t' * G * obj.Dm_inv[2*t-1:2*t,:]
    end
end

function update_mesh(mesh::Mesh, obj::CGTriObject)
    for i in 1:obj.N
        mesh.vertices[i].x = obj.x_node[i,:]
    end
end
