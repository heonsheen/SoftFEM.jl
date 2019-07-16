include("ElasticObject.jl")
using SparseArrays
using SharedArrays
using Distributed
using ExpmV

# operator splitting with IM and EX, parallel solve elasticity
# TODO: implement ERE
function dg_mixed_integrator(u::Vector{Float64},
                            obj::ElasticObject,
                            dt::Float64,
                            fixed::Vector{Bool},
                            g::Vector{Float64}, # gravity / external force
                            els_int::String = "IM", # elasticity integrator
                            int_int::String = "EX") # interface integrator
    free_ind = .!fixed;

    n = obj.N * obj.dim
    n_fixed = convert(Int64, sum(fixed) / obj.dim)
    n_free = convert(Int64, sum(free_ind) / obj.dim)
    q = u[1:n]
    v = u[n+1:end]
    u_new = deepcopy(u)
    q_new = u_new[1:n]
    v_new = u_new[n+1:end]

    f_cnt = 0
    f_inds = []
    for i = 1:obj.NT
        inds = [obj.dim*(j-1)+k for k in 1:obj.dim, j in obj.ec[i,:]]
        n_ft = convert(Int64, sum(free_ind[inds]))
        push!(f_inds, f_cnt+1:f_cnt+n_ft)
        f_cnt += n_ft
    end

    M = obj.M
    M = M[free_ind,free_ind]

    obj.x = obj.X + q_new + v_new * dt
    update_pos(obj, obj.x - obj.X)

    K_int = compute_interface_stiffness_matrix(obj)
    K_int = K_int[free_ind, free_ind]

    f_int = compute_interface_force(obj)
    f_int = f_int[free_ind]

    B_int = obj.mat.alpha * M + obj.mat.beta * K_int
    f_2 = f_int - B_int*(v_new[free_ind])

    Dv_int = zeros(Float64,n)
    if int_int == "EX"
        Dv_int[free_ind] = dt * (M \ f_int)
    elseif int_int == "ERE"
        J = [spzeros(obj.dim*n_free,obj.dim*n_free) sparse(I,obj.dim*n_free,obj.dim*n_free);
            -M\K_int -M\B_int]
        du = [v[free_ind]; M\f_2]
        c_n = du - J*[q[free_ind]; v[free_ind]]
        nJ = 2*obj.dim*n_free
        A = [J c_n; spzeros(1,nJ+1)]
        u_tilde = [q[free_ind]; v[free_ind]; 1]

        #u_p = [sparse(I,nJ,nJ) spzeros(nJ,1)] * exp(dt*A) * u_tilde
        I0 = [sparse(I,nJ,nJ) spzeros(nJ,1)] 
        u_p = I0 * real(expmv(dt, A, u_tilde))
        Dv_int[free_ind] = u_p[obj.dim*n_free+1:end] - v[free_ind]
    end

    max_iters = 20
    iter = 0

    if els_int == "IM"
        # Newton steps to solve the system
        # TODO: Implement line search
        while true
            obj.x = obj.X + q_new + v_new * dt + Dv_int*dt
            update_pos(obj, obj.x - obj.X)

            K_els = compute_elastic_stiffness_matrix(obj)
            f_els = compute_elastic_force(obj)

            # stiffness matrix
            K_els = K_els[free_ind,free_ind]

            # damping
            B = obj.mat.alpha * M + obj.mat.beta * K_els

            # force (RHS)
            f_els = f_els[free_ind]
            f_ext = M * g[free_ind]
            f_1 = f_els + f_ext - B*(v_new[free_ind])
            
            # parallel block diag
            Dv = zeros(Float64, n_free * obj.dim)
            Dvs = convert(SharedArray, Dv)
            @sync @distributed for i in 1:obj.NT
                inds = f_inds[i]
                LHS = I + dt*dt*(M[inds,inds]\K_els[inds,inds]) + 
                            dt*(M[inds,inds]\(B[inds,inds]))
                RHS = v_new[free_ind][inds] - v[free_ind][inds] - 
                        Dv_int[free_ind][inds] - dt*(M[inds,inds]\f_1[inds])
                Dvs[inds] = -LHS \ RHS
            end
            Dv = convert(Array, Dvs)
            
            #=
            LHS = sparse(I,obj.dim*n_free,obj.dim*n_free) + dt*dt*(M\K_els) - dt*(M\B)
            RHS = v_new[free_ind] - v[free_ind] - Dv_int[free_ind] - dt*(M\f_1)
            Dv = -LHS \ RHS
            =#
            v_new[free_ind] += Dv
            
            residual = (v_new[free_ind] - v[free_ind] - Dv_int[free_ind] - dt * (M\f_1))' * 
                        (v_new[free_ind] - v[free_ind] - Dv_int[free_ind] - dt * (M\f_1))
            iter = iter + 1
            
            u_new[1:n] = q + dt * (v_new)
            u_new[n+1:end] = v_new

            #println("Dv^T Dv = ", Dv'*Dv)
            #println("residual = ", residual)

            ((Dv'*Dv <= 1e-12) || (residual <= 1e-12) || iter >= max_iters) && break
        end
        
        println("# iters = ", iter)
    end

    obj.x = obj.X + q_new + v_new * dt
    obj.v = v_new
    update_pos(obj, obj.x - obj.X)

    u_new
end