include("ElasticObject.jl")
include("Geometry.jl")
include("Material.jl")
include("Integrals.jl")

# Discontinous Galerkin method with Triangular Finite Elements
mutable struct DGTriObject <: ElasticObject
### Attributes 
    N::Int64 # number of nodes
    NT::Int64 # number of elements
    dim::Int64 # dimension of the object

    x_node::Matrix{Float64} # nodal position in world coordinates (N x dim)
    X_node::Matrix{Float64} # nodal position in material coordinates (N x dim)
    ec::Matrix{Int64} # element connectivity (NT x {3,4,6,...})

    N_CG::Int64 # number of CG nodes
    ec_CG::Matrix{Int64} # CG element connectivity
    DG_map::Vector{Int64} # Maps DG element indices to CG element indices

    NE::Int64 # number of interface edges
    int_minus_edges::Array{Tuple{Int64,Int64}} # DG Interface edge list on the minus side
    int_plus_edges::Array{Tuple{Int64,Int64}} # DG Interface edge list on the plus side
    interface_elem::Array{Tuple{Int64,Int64}} # DG Interface adjacent face indices per edge

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

    L::Vector{Float64} # reference lengths of each interface edges (NE x 1)
    nor::Vector{Float64} # unit length outward normals of interface edges ((dim * NE) x 1)

    b::Vector{Float64} # Translation components need for DG ((dim * NT) x 1)

    T::Matrix{Float64} # mapping vectorized nodal position in a tri to 
                       #   its vectorized deformation gradient (4NT by 6)
                       #   definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)

    M::SparseMatrixCSC{Float64,Int64} # Mass matrix
    K_prev::SparseMatrixCSC{Float64,Int64} # Stiffness matrix of previous timestep
    f_prev::SparseMatrixCSC{Float64,Int64} # force vector of previous timestep
    K0::SparseMatrixCSC{Float64,Int64} # 

    mat::Material # elastic material description
    
### Constructor
    function DGTriObject(
        mesh::Mesh,
        mat::Material)
        NT = mesh.n_faces
        N = NT * 3
        N_CG = mesh.n_vertices
        dim = mesh.dim
        ec_CG = mesh.ec

        ec = zeros(Int64, NT, 3)
        DG_map = zeros(Int64, NT*3)
        for t in 1:NT
            ec[t,:] = [3*(t-1)+1, 3*(t-1)+2, 3*t]
            for j in 1:3
                DG_map[3*(t-1)+j] = mesh.ec[t,j]
            end
        end

        x_node = [mesh.vertices[DG_map[i]].x[j] for i in 1:N, j = 1:dim]
        X_node = x_node

        x = vec(reshape(x_node', (dim*N, 1)))
        X = vec(reshape(x_node', (dim*N, 1)))
        v = zeros(dim*N)

        interface_elem = []
        int_minus_edges = []
        int_plus_edges = []
        L = []
        nor = []
        for t in 1:NT, i in 1:3
            te = mesh.e2e[t, i]
            # edge is not on boundary && edge not already in elem list
            if te[2] != 0 && !((te[1],t) in interface_elem) 
                interface_elem = [interface_elem; (t,te[1])]
                he_m = mesh.half_edges[te[1]][te[2]]
                he_p = mesh.half_edges[t][i]
                @assert he_m.origin == he_p.dest && he_m.dest == he_p.origin "Half edge pairing incorrect"
                int_minus_edges = [int_minus_edges; 
                                    (ec[te[1],te[2]], ec[te[1],mod(te[2],3)+1])]
                int_plus_edges = [int_plus_edges; 
                                    (ec[t,i], ec[t,mod(i,3)+1])]
                X0 = mesh.vertices[he_p.origin].x
                X1 = mesh.vertices[he_p.dest].x
                L = [L; norm(X0-X1)]

                e_p = X1 - X0
                n_p = [e_p[2], -e_p[1]]
                n_p = normalize(n_p)
                nor = [nor; n_p]
            end
        end
        NE = size(interface_elem,1)
        
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
        b = zeros(dim*NT)
        
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
        
        new(N, NT, dim, x_node, X_node, ec, 
            N_CG, ec_CG, DG_map, NE, int_minus_edges, int_plus_edges, interface_elem,
            x, X, v, Ds, Dm, Dm_inv, F, F_inv, dF, W, L, nor, b, T, 
            M, K_prev, f_prev, K0, mat)
    end
end

function compute_elastic_stiffness_matrix(obj::DGTriObject)
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

function compute_elastic_force(obj::DGTriObject)
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

function compute_interface_stiffness_matrix(obj::DGTriObject)
    # BZ implementation
    # TODO: move BZ, IP specification to different struct/type
    I = zeros(Int64,9*4*obj.NT)
    J = zeros(Int64,9*4*obj.NT)
    V = zeros(Float64,9*4*obj.NT)
    nnz = 0

    for e in 1:obj.NE
        ef = obj.interface_elem[e]
        fi_m = ef[1] # face index minus
        fi_p = ef[2] # face index plus

        F_m = obj.F[2*(fi_m-1)+1:2*fi_m,:]
        F_p = obj.F[2*(fi_p-1)+1:2*fi_p,:]

        b_m = obj.b[2*(fi_m-1)+1:2*fi_m]
        b_p = obj.b[2*(fi_p-1)+1:2*fi_p]

        Dmi_m = obj.Dm_inv[2*(fi_m-1)+1:2*fi_m,:]
        Dmi_p = obj.Dm_inv[2*(fi_p-1)+1:2*fi_p,:]

        len = obj.L[e]

        v0_m = obj.int_minus_edges[e][1]
        v1_m = obj.int_minus_edges[e][2]
        v0_p = obj.int_plus_edges[e][1]
        v1_p = obj.int_plus_edges[e][2]

        P0 = obj.X_node[v0_p,:]
        P1 = obj.X_node[v1_p,:]

        n_p = obj.nor[2*(e-1)+1:2*e]

        I2 = Matrix{Float64}(I, 2, 2)
        G = [1 0; 0 1; -1 -1]
        bary = [1 0 0]
        eta_t = obj.mat.eta * len * (1.0/obj.W[fi_m] + 1.0/obj.W[fi_p]);

        for i = 1:3, j = 1:2, k = 1:3, l = 1:2
            # minus side element
            dFj_m = I2[:,j] * G[i,:]' * Dmi_m
            dFl_m = I2[:,l] * G[k,:]' * Dmi_m
            dbj_m = I2[:,j] * bary[i] - 
                I2[:,j] * G[i,:]' * Dmi_m * obj.X_node[obj.ec[fi_m,1],:]
            dbl_m = I2[:,l] * bary[k] - 
                I2[:,l] * G[k,:]' * Dmi_m * obj.X_node[obj.ec[fi_m,1],:]

            # plus side element
            dFj_p = I2[:,j] * G[i,:]' * Dmi_p
            dFl_p = I2[:,l] * G[k,:]' * Dmi_p
            dbj_p = I2[:,j] * bary[i] - 
                I2[:,j] * G[i,:]' * Dmi_p * obj.X_node[obj.ec[fi_p,1],:]
            dbl_p = I2[:,l] * bary[k] - 
                I2[:,l] * G[k,:]' * Dmi_p * obj.X_node[obj.ec[fi_p,1],:]    
            
            # d/dx-^l, d/dx-^j E    
            A_mm = dFj_m' * dFl_m + dFl_m' * dFj_m
            B_mm = 2 * dbl_m' * dFj_m + 2 * dbj_m' * dFl_m
            c_mm = dbj_m' * dbl_m + dbl_m' * dbj_m

            ddE_mm = quadratic_line_int(A_mm, P0, P1) + 
                    linear_line_int(Matrix{Float64}(B_mm), P0, P1) + 
                    const_line_int(c_mm, P0, P1)
            push!(I, 2*(obj.ec[fi_m,k]-1)+l)
            push!(J, 2*(obj.ec[fi_m,i]-1)+j)
            push!(V, 0.5 * eta_t * ddE_mm)

            # d/dx+^l, d/dx-^j E  
            A_pm = -dFj_m' * dFl_p - dFl_p' * dFj_m
            B_pm = 2 * -dbl_p' * dFj_m + 2 * dbj_m' * -dFl_p
            c_pm = dbj_m' * -dbl_p + -dbl_p' * dbj_m

            ddE_pm = quadratic_line_int(A_pm, P0, P1) + 
                    linear_line_int(Matrix{Float64}(B_pm), P0, P1) + 
                    const_line_int(c_pm, P0, P1)
            push!(I, 2*(obj.ec[fi_p,k]-1)+l)
            push!(J, 2*(obj.ec[fi_m,i]-1)+j)
            push!(V, 0.5 * eta_t * ddE_pm)
        
            # d/dx-^l, d/dx+^j E  
            A_mp = -dFj_p' * dFl_m - dFl_m' * dFj_p
            B_mp = 2 * dbl_m' * -dFj_p + 2 * -dbj_p' * dFl_m
            c_mp = -dbj_p' * dbl_m + dbl_m' * -dbj_p

            ddE_mp = quadratic_line_int(A_mp, P0, P1) + 
                    linear_line_int(Matrix{Float64}(B_mp), P0, P1) + 
                    const_line_int(c_mp, P0, P1)
            push!(I, 2*(obj.ec[fi_m,k]-1)+l)
            push!(J, 2*(obj.ec[fi_p,i]-1)+j)
            push!(V, 0.5 * eta_t * ddE_mp)

            # d/dx+^l, d/dx+^j E  
            A_pp = dFj_p' * dFl_p + dFl_p' * dFj_p
            B_pp = 2 * -dbl_p' * -dFj_p + 2 * -dbj_p' * -dFl_p
            c_pp = -dbj_p' * -dbl_p + -dbl_p' * -dbj_p

            ddE_pp = quadratic_line_int(A_pp, P0, P1) + 
                    linear_line_int(Matrix{Float64}(B_pp), P0, P1) + 
                    const_line_int(c_pp, P0, P1)
            push!(I, 2*(obj.ec[fi_p,k]-1)+l)
            push!(J, 2*(obj.ec[fi_p,i]-1)+j)
            push!(V, 0.5 * eta_t * ddE_pp)
        end
    end

    sparse(I, J, V, 2*obj.N, 2*obj.N)
end

function compute_interface_force(obj::DGTriObject)
    # BZ implementation
    # TODO: move BZ, IP specification to different struct/type
    f = zeros(Float64, obj.dim * obj.N)

    for e in 1:obj.NE
        ef = obj.interface_elem[e]
        fi_m = ef[1] # face index minus
        fi_p = ef[2] # face index plus

        F_m = obj.F[2*(fi_m-1)+1:2*fi_m,:]
        F_p = obj.F[2*(fi_p-1)+1:2*fi_p,:]

        b_m = obj.b[2*(fi_m-1)+1:2*fi_m]
        b_p = obj.b[2*(fi_p-1)+1:2*fi_p]

        Dmi_m = obj.Dm_inv[2*(fi_m-1)+1:2*fi_m,:]
        Dmi_p = obj.Dm_inv[2*(fi_p-1)+1:2*fi_p,:]

        len = obj.L[e]

        v0_m = obj.int_minus_edges[e][1]
        v1_m = obj.int_minus_edges[e][2]
        v0_p = obj.int_plus_edges[e][1]
        v1_p = obj.int_plus_edges[e][2]

        P0 = obj.X_node[v0_p,:]
        P1 = obj.X_node[v1_p,:]

        n_p = obj.nor[2*(e-1)+1:2*e]

        I2 = Matrix{Float64}(I, 2, 2)
        G = [1 0; 0 1; -1 -1]
        bary = [1 0 0]
        eta_t = obj.mat.eta * len * (1.0/obj.W[fi_m] + 1.0/obj.W[fi_p]);

        for i = 1:3, j = 1:2
            # minus side element
            dF_m = I2[:,j] * G[i,:]' * Dmi_m
            db_m = I2[:,j] * bary[i] - 
                I2[:,j] * G[i,:]' * Dmi_m * obj.X_node[obj.ec[fi_m,1],:]
            
            A_m = dF_m' * F_m + F_m' * dF_m - dF_m' * F_p - F_p' * dF_m
            B_m = 2 * (b_m - b_p)' * dF_m + 2 * db_m' * (F_m - F_p)
            c_m = db_m' * (b_m - b_p) + (b_m - b_p)' * db_m

            dE_m = quadratic_line_int(A_m, P0, P1) + 
                    linear_line_int(Matrix{Float64}(B_m), P0, P1) + 
                    const_line_int(c_m, P0, P1)
            f[2*(obj.ec[fi_m,i]-1)+j] += -0.5 * eta_t * dE_m

            # plus side element
            dF_p = I2[:,j] * G[i,:]' * Dmi_p
            db_p = I2[:,j] * bary[i] - 
                I2[:,j] * G[i,:]' * Dmi_p * obj.X_node[obj.ec[fi_p,1],:]
            
            A_p = dF_p' * F_p + F_p' * dF_p - dF_p' * F_m - F_m' * dF_p
            B_p = 2 * (b_m - b_p)' * -dF_p + 2 * -db_p' * (F_m - F_p)
            c_p = -db_p' * (b_m - b_p) + (b_m - b_p)' * -db_p

            dE_p = quadratic_line_int(A_p, P0, P1) + 
                    linear_line_int(Matrix{Float64}(B_p), P0, P1) + 
                    const_line_int(c_p, P0, P1)
            f[2*(obj.ec[fi_p,i]-1)+j] += -0.5 * eta_t * dE_p
        end
    end

    f
end

function compute_force_differential(obj::DGTriObject)
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

function update_pos(obj::DGTriObject, dx::Vector{Float64})
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
        X1 = obj.X_node[obj.ec[t,1],:]
        x1 = obj.x_node[obj.ec[t,1],:]
        obj.b[2*t-1:2*t] = x1 - obj.F[2*t-1:2*t,:] * X1
    end
end
#=
function update_mesh(mesh::Mesh, obj::DGTriObject)
    for i in 1:obj.N
        mesh.vertices[i].x = obj.x_node[i,:]
    end
end
=#
function get_DG_mesh(obj::DGTriObject)
    vertices = Array{Vertex}(undef, obj.N)
    for i in 1:obj.N
        vertices[i] = Vertex(obj.x_node[i,:])
    end
    mesh = Mesh(vertices, obj.ec)
    mesh
end
