include("ElasticObject.jl")
using SparseArrays

function forward_euler(u::Vector{Float64},
                        obj::ElasticObject,
                        dt::Float64,
                        fixed::Vector{Bool},
                        g::Vector{Float64}) # gravity / external force
    free_ind = .!fixed;

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

    # force (RHS)
    f_el = f_el[free_ind]
    f_ext = M * g[free_ind]
    f = f_el + f_ext
            
    v_new[free_ind] += dt * (M\f)
    
    u_new[1:n] = q + dt * (v_new)
    u_new[n+1:end] = v_new
    
    obj.x = obj.X + q_new + v_new * dt
    obj.v = v_new
    update_pos(obj, obj.x - obj.X)

    u_new
end