include("ElasticObject.jl")
using SparseArrays
using ExpmV

function ERE(u::Vector{Float64},
                obj::ElasticObject,
                dt::Float64,
                fixed::Vector{Bool},
                g::Vector{Float64}) # gravity / external force
    free_ind = .!fixed

    n = obj.N * obj.dim
    n_fixed = convert(Int64, sum(fixed) / obj.dim)
    n_free = convert(Int64, sum(free_ind) / obj.dim)
    q = u[1:n]
    v = u[n+1:end]
    u_new = deepcopy(u)
    q_new = u_new[1:n]
    v_new = u_new[n+1:end]

    M = obj.M
    M = M[free_ind,free_ind]

    obj.x = obj.X + q_new + v_new * dt
    update_pos(obj, obj.x - obj.X)

    K = compute_total_stiffness_matrix(obj)
 
    f_el = compute_total_force(obj)

    # stiffness matrix
    K = K[free_ind,free_ind]

    K0 = obj.K0[free_ind, free_ind]

    # damping
    B = obj.mat.alpha * M + obj.mat.beta * K0

    # force (RHS)
    f_el = f_el[free_ind]
    f_ext = M * g[free_ind]
    f = f_el + f_ext - B*(v_new[free_ind])

    J = [spzeros(obj.dim*n_free,obj.dim*n_free) sparse(I,obj.dim*n_free,obj.dim*n_free);
        -M\K -M\B]
    du = [v[free_ind]; M\f]
    c_n = du - J*[q[free_ind]; v[free_ind]]
    nJ = 2*obj.dim*n_free
    A = [J c_n; spzeros(1,nJ+1)]
    u_tilde = [q[free_ind]; v[free_ind]; 1]

    #u_p = [sparse(I,nJ,nJ) spzeros(nJ,1)] * exp(dt*A) * u_tilde
    I0 = [sparse(I,nJ,nJ) spzeros(nJ,1)] 
    u_p = I0 * expmv(dt, A, u_tilde)
    #=
    v_new[free_ind] += Dv
    u_new[1:n] = q + dt * (v_new)
    u_new[n+1:end] = v_new
    =#
    q_new[free_ind] = u_p[1:obj.dim*n_free]
    v_new[free_ind] = u_p[obj.dim*n_free+1:end]
    u_new = [q_new; v_new]

    obj.x = obj.X + q_new + v_new * dt
    obj.v = v_new
    update_pos(obj, obj.x - obj.X)

    u_new
end