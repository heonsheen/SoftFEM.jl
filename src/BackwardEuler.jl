include("ElasticObject.jl")
using SparseArrays

function backward_euler(u::Vector{Float64},
                        obj::ElasticObject,
                        dt::Float64,
                        fixed::Vector{Bool},
                        g::Vector{Float64}) # gravity / external force
    free_ind = .!fixed;

    n = obj.N * obj.dim
    n_fixed = convert(Int64, sum(fixed) / obj.dim)
    n_free = convert(Int64, sum(free_ind) / obj.dim)
    dx = u[1:n]
    v = u[n+1:end]
    u_new = deepcopy(u)
    dx_new = u_new[1:n]
    v_new = u_new[n+1:end]

    M = obj.M
    M = M[free_ind,free_ind]

    max_iters = 20
    iter = 0

    while true
        obj.x = obj.X + dx_new + v_new * dt
        update_pos(obj, obj.x - obj.X)

        K = compute_stiffness_matrix(obj)
        f_el = compute_elastic_force(obj)

        K = K[free_ind,free_ind]
        f_el = f_el[free_ind]
        f_ext = M * g[free_ind]
        f = f_el + f_ext

        B = 0 * M - 0 * K # TODO: make this damping

        Dv = -(sparse(I,obj.dim*n_free,obj.dim*n_free) + dt*dt*(M\K) - dt*(M\B)) \ 
                (v_new[free_ind] - v[free_ind] - dt*(M\f))
        v_new[free_ind] += Dv

        residual = (v_new[free_ind] - v[free_ind] - dt * (M\f))' * 
                    (v_new[free_ind] - v[free_ind] - dt * (M\f))
        iter = iter + 1
        
        u_new[1:n] = dx + dt * (v_new)
        u_new[n+1:end] = v_new

        # println("Dv^T Dv = ", Dv'*Dv)
        # println("residual = ", residual)
        ((Dv'*Dv <= 1e-12) || (residual <= 1e-12) || iter >= max_iters) && break
    end
    
    obj.x = obj.X + dx_new + v_new * dt
    obj.v = v_new
    update_pos(obj, obj.x - obj.X)

    u_new
end