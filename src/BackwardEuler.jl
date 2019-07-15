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
    q = u[1:n]
    v = u[n+1:end]
    u_new = deepcopy(u)
    q_new = u_new[1:n]
    v_new = u_new[n+1:end]

    M = obj.M
    M = M[free_ind,free_ind]

    max_iters = 20
    iter = 0

    # Newton steps to solve the system
    while true
        obj.x = obj.X + q_new + v_new * dt
        update_pos(obj, obj.x - obj.X)

        K = compute_total_stiffness_matrix(obj)
        f_el = compute_total_force(obj)

        # stiffness matrix
        K = K[free_ind,free_ind]

        # damping
        B = obj.mat.alpha * M - obj.mat.beta * K

        # force (RHS)
        f_el = f_el[free_ind]
        f_ext = M * g[free_ind]
        f = f_el + f_ext + B*(v_new[free_ind])

        LHS = sparse(I,obj.dim*n_free,obj.dim*n_free) + dt*dt*(M\K) - dt*(M\B)
        RHS = v_new[free_ind] - v[free_ind] - dt*(M\f)
        Dv = -LHS \ RHS
                
        v_new[free_ind] += Dv

        residual = (v_new[free_ind] - v[free_ind] - dt * (M\f))' * 
                    (v_new[free_ind] - v[free_ind] - dt * (M\f))
        iter = iter + 1
        
        u_new[1:n] = q + dt * (v_new)
        u_new[n+1:end] = v_new

        println("Dv^T Dv = ", Dv'*Dv)
        println("residual = ", residual)
        ((Dv'*Dv <= 1e-12) || (residual <= 1e-12) || iter >= max_iters) && break
    end

    println("# iters = ", iter)
    
    obj.x = obj.X + q_new + v_new * dt
    obj.v = v_new
    update_pos(obj, obj.x - obj.X)

    u_new
end