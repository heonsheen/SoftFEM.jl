# A simple tet mesh structure

mutable struct Node3d{T<:Real}
    # Reference Positions
    X::Tuple{T,T,T}
    # Current Positions
    x::Tuple{T,T,T}
    # Constructor
    function Node3d{T}(X::Tuple{T,T,T}, x::Tuple{T,T,T}) where T
        new(X, x)
    end

    function Node3d{T}{X::Tuple{T,T,T}) where T
        new(X, X)
    end
end

mutable struct TetElem{T<:Real}
    # Element Topology
    topo::Tuple{Int64,Int64,Int64,Int64}
    # Element Nodes
    nodes::Tuple{Node3d{T},Node3d{T},Node3d{T},Node3d{T}}
    # Constructor
    function TetElem(
        topo::Tuple{Int64,Int64,Int64,Int64},
        nodes::Tuple{Node3d{T},Node3d{T},Node3d{T},Node3d{T}} ) where T
        new(topo, nodes)
    end
end

mutable struct Tetmesh{T<:Real}
    n_nodes::Int64
    n_elems::Int64
    nodes::Array{Node3d{T}}
    elems::Array{TetElem{T}}
    internal_nodes_topo::Array{Int64}
    boundary_nodes_topo::Array{Int64}
    adjacency_list::Array{Int64,2}

    function Tetmesh(
        n_nodes::Int64,
        n_elems::Int64,
        nodes::Array{Node3d{T}},
        elems::Array{TetElem{T}} ) where T
        @assert(size(nodes) == n_nodes, 'Nodes size not consistent' )
        @assert(size(elems) == n_elems, 'Elements size not consistent' )
        
        for i = 1 : n_elems
            T_i = elems[i]
            al_i = []
            for j = 1 : n_elems
                if i != j
                    T_j = elems[j]
                    for k = 1 : 4
                        tri_topo = T_i.topo[1:end .!= k]
                        if (tri_topo[1] in T_j) && (tri_topo[2] in T_j) && (tri_topo[3] in T_j)
                            push!(al_i, j)
                            for l = 1 : 3
                                if !(tri_topo[l] in internal_nodes_topo)
                                    push!(internal_nodes_topo, tri_topo[l])
                                end
                            end
                        end
                    end
                end
            end
            push!(adjacency_list, al_i)
        end

        boundary_nodes_topo = setdiff(nodes, internal_nodes_topo)
    end
    
end

# A simple tet mesh structure

mutable struct Node2d{T<:Real}
    # Reference Positions
    X::Tuple{T,T}
    # Current Positions
    x::Tuple{T,T}
    # Constructor
    function Node3d{T}(X::Tuple{T,T}, x::Tuple{T,T}) where T
        new(X, x)
    end

    function Node3d{T}{X::Tuple{T,T}) where T
        new(X, X)
    end
end

mutable struct TriElem{T<:Real}
    # Element Topology
    topo::Tuple{Int64,Int64,Int64}
    # Element Nodes
    nodes::Tuple{Node2d{T},Node2d{T},Node2d{T}}
    # Constructor
    function TriElem(
        topo::Tuple{Int64,Int64,Int64},
        nodes::Tuple{Node2d{T},Node2d{T},Node2d{T}} ) where T
        new(topo, nodes)
    end
end

mutable struct Trimesh{T<:Real}
    n_nodes::Int64
    n_elems::Int64
    nodes::Array{Node2d{T}}
    elems::Array{TriElem{T}}
    internal_nodes_topo::Array{Int64}
    boundary_nodes_topo::Array{Int64}
    adjacency_list::Array{Int64,2}

    function Trimesh(
        n_nodes::Int64,
        n_elems::Int64,
        nodes::Array{Node2d{T}},
        elems::Array{TriElem{T}} ) where T
        @assert(size(nodes) == n_nodes, 'Nodes size not consistent' )
        @assert(size(elems) == n_elems, 'Elements size not consistent' )
        
        for i = 1 : n_elems
            T_i = elems[i]
            al_i = []
            for j = 1 : n_elems
                if i != j
                    T_j = elems[j]
                    for k = 1 : 3
                        tri_topo = T_i.topo[1:end .!= k]
                        if (tri_topo[1] in T_j) && (tri_topo[2] in T_j)
                            push!(al_i, j)
                            for l = 1 : 2
                                if !(tri_topo[l] in internal_nodes_topo)
                                    push!(internal_nodes_topo, tri_topo[l])
                                end
                            end
                        end
                    end
                end
            end
            push!(adjacency_list, al_i)
        end

        boundary_nodes_topo = setdiff(nodes, internal_nodes_topo)
    end
    
end
