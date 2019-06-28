include("ElasticObject.jl")

function backward_euler(u::Vector{Float64},
                        obj::ElasticObject,
                        dt::Float64,
                        fixed::Vector{Bool},
                        g::Vector{Float64}) # gravity / external force
    free_ind = .!fixed;

    N = obj.N * obj.dim
    dx = u[1:N]
    v = u[N+1:end]
    u_new = deepcopy(u)
    dx_new = u_new[1:N]
    v_new = u_new[N+1:end]

    obj.x = obj.X + dx + v * dt
    update_pos(obj, dx)

    K = compute_stiffness_matrix(obj)
    f_el = compute_elastic_force(obj)
    M = obj.M

    K = K[free_ind,free_ind]
    M = M[free_ind,free_ind]
    f_el = f_el[free_ind]
    f_ext = M * g[free_ind]
    f = f_el + f_ext

    n_fixed = sum(fixed) / obj.dim
    v_new[free_ind] += -(sparse(I,obj.dim*(N-n_fixed)) + dt*dt*(M\K))

    
end