include("ElasticObject.jl")
using SparseArrays

function backward_euler(u::Vector{Float64},
                        obj::ElasticObject,
                        dt::Float64,
                        fixed::Vector{Bool},
                        g::Vector{Float64},# gravity / external force
                        line_search::Bool = false) 
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
    K = spzeros(n,n)
    f = zeros(n) 
    
    K0 = obj.K0[free_ind,free_ind]
    B = obj.mat.alpha * M + obj.mat.beta * K0
    
    println("norm(B) = ", norm(B))
    
    f_ext = M * g[free_ind]

    max_iters = 40
    iter = 0

    # Newton steps to solve the system
    # TODO: Implement line search
    while true
        obj.x = obj.X + q_new + v_new * dt
        update_pos(obj, obj.x - obj.X)

        K = compute_total_stiffness_matrix(obj)
        f_el = compute_total_force(obj)

        # stiffness matrix
        K = K[free_ind,free_ind]

        # damping
        #B = obj.mat.alpha * M + obj.mat.beta * K

        # force (RHS)
        f_el = f_el[free_ind]
        f = f_el + f_ext - B*(v_new[free_ind])

        LHS = sparse(I,obj.dim*n_free,obj.dim*n_free) + dt*dt*(M\K) + dt*(M\B)
        RHS = v_new[free_ind] - v[free_ind] - dt*(M\f)
        Dv = -LHS \ RHS
                
        if !line_search
            v_new[free_ind] += Dv

            residual = (v_new[free_ind] - v[free_ind] - dt * (M\f))' * 
                        (v_new[free_ind] - v[free_ind] - dt * (M\f))
            iter = iter + 1
            
            u_new[1:n] = q + dt * (v_new)
            u_new[n+1:end] = v_new
        else
            residual_old = 0.5 *(dt * (M\f))' * (dt * (M\f))
            
            step_size = 1.0
            v_new_temp = v_new
            c1 = 1e-4
            c2 = 0.9
            
            v_new_temp[free_ind] = v_new[free_ind] + step_size*Dv
            q_new_temp = q + dt * (v_new_temp)
            update_pos(obj, q_new_temp)
            K = compute_total_stiffness_matrix(obj)
            K = K[free_ind, free_ind]
            #B = obj.mat.alpha * M + obj.mat.beta * K
            f_els = compute_total_force(obj)
            f_els = f_els[free_ind]
            f = f_els + f_ext - B*v_new_temp[free_ind]
            
            # merit function
            residual = 0.5*(v_new_temp[free_ind] - v[free_ind] - dt * (M\f))' * 
                            (v_new_temp[free_ind] - v[free_ind] - dt * (M\f))
            ls_it = 0
            
            # wolfe conditions
            while (residual > residual_old - c1 * step_size * residual * 2)
                residual_old = residual
                ls_it = ls_it + 1
                if ls_it == 40
                    error("line search more than 40 iterations.")
                end
                step_size = 0.9 * step_size
                v_new_temp[free_ind] = v_new[free_ind] + step_size*Dv
                q_new_temp = q + dt * (v_new_temp)
                update_pos(obj, q_new_temp)
                K = compute_total_stiffness_matrix(obj)
                K = K[free_ind, free_ind]
                #B = obj.mat.alpha * M + obj.mat.beta * K
                f_els = compute_total_force(obj)
                f_els = f_els[free_ind]
                f = f_els + f_ext - B*v_new_temp[free_ind]
                
                # merit function
                residual = 0.5*(v_new_temp[free_ind] - v[free_ind] - dt * (M\f))' * 
                                (v_new_temp[free_ind] - v[free_ind] - dt * (M\f))
                println("line search residual = ", residual)
            end

            v_new[free_ind] = v_new_temp[free_ind]

            residual = (v_new[free_ind] - v[free_ind] - dt * (M\f))' * 
                        (v_new[free_ind] - v[free_ind] - dt * (M\f))
            iter = iter + 1
            
            u_new[1:n] = q + dt * (v_new)
            u_new[n+1:end] = v_new
            println("line search # iters = ", ls_it)
        end

        println("Dv^T Dv = ", Dv'*Dv)
        println("residual = ", residual)
        ((Dv'*Dv <= 1e-12) || (residual <= 1e-12) || iter >= max_iters) && break
    end

    println("# iters = ", iter)
    
    obj.x = obj.X + q_new + v_new * dt
    obj.v = v_new
    update_pos(obj, obj.x - obj.X)

    obj.K_prev = K
    obj.f_prev = f

    u_new
end