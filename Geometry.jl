# Simple definitions for basic geometry structure
# Following Alumbaugh et al. 2005 (Compact Array-Based Mesh Data Structures)

using StaticArrays
using LinearAlgebra

# N-dimensional Vertex data
mutable struct Vertex
    ### Attributes
    # Positions of the vertex
    x::Array{Float64}

    ### Constructors
    function Vertex(x::Array{Float64})
        N = size(x,1)
        @assert N >= 2 "Dimension cannot be less than 2"
    	@assert N <= 3 "Dimension cannot be larger than 3"
    	
        new(x)
    end
end

# Half-edge connectivity data
mutable struct HalfEdge
    ### Attributes
    # Index of the vertex the half-edge originates from
    origin::Int64
    # Index of the vertex the half-edge points to
    dest::Int64

    ### Constructors
    function HalfEdge(
	origin::Int64,
        dest::Int64)
	if origin < 0 || dest < 0
	    error("Invalid indexing.")
	end
	new(origin, dest)
    end
end

# Half-Face connectivity data for volumetric mesh topology
mutable struct HalfFace
    ### Attributes
    # Indices of the vertices belonging to the half-face
    vertices::Array{Int64}
    # Local index of the anchor vertex
    anchor::Int64

    ### Constructor
    function HalfFace(vertices::Array{Int64}, anchor::Int64)
	new(vertices, anchor)
    end
end

# Uniform surface mesh structure with dimensionality N
# TODO: Support non-uniform mesh
mutable struct Mesh
    ### Attributes
    # Number of vertices in the mesh
    n_vertices::Int64
    # Number of faces in the mesh
    n_faces::Int64
    # Number of boundary half-edges
    n_boundary::Int64
    # Array of vertices in the mesh
    vertices::Array{Vertex}
    # 2D Array of half-edges in the mesh indexed by [face id, local id])
    half_edges::Array{Array{HalfEdge}}
    # Array of boundary half-edges
    # Element Connectivity Array
    ec::Array{Int64,2}
    # Map of vertex to half-edge originating from it
    #   & Map of boundary vertex to boundary half-edge
    v2e::Array{Tuple{Int64,Int64}}
    # Map of internal half-edge to its twin half-edge
    e2e::Array{Tuple{Int64,Int64},2}
    # Map of boundary half-edge to its twin internal half-edge
    b2e::Array{Tuple{Int64,Int64}}
    
    ### Constructor & check mesh sanity
    function Mesh(
        vertices::Array{Vertex},
        ec::Array{Int64,2} )
    	# populate number of elements
        n_vertices = size(vertices,1)
        n_faces = size(ec,1)
        elem_dim = size(ec,2)
        
        # Populate half-edge array
        vert_he_map = Tuple{Int64,Int64}[(0,0) for i in 1:n_vertices, j in 1:n_vertices]
        half_edges = Vector{Vector{HalfEdge}}()
        for fid in 1:n_faces
            face_half_edges = Array{HalfEdge}(undef, elem_dim)
            for lid in 1:elem_dim
                he = HalfEdge(ec[fid, lid], ec[fid, mod(lid, elem_dim)+1])
                face_half_edges[lid] = he
                vert_he_map[ec[fid,lid], ec[fid,mod(lid, elem_dim)+1]] = (fid, lid)
            end
            push!(half_edges, face_half_edges)
        end

        # Populate half-edge - twin half-edge map
        e2e = Tuple{Int64,Int64}[(0,0) for i in 1:n_faces, j in 1:elem_dim]
        v2e = Tuple{Int64,Int64}[(0,0) for i in 1:n_vertices]
        b2e = Tuple{Int64,Int64}[]
        n_boundary = 0 # boundary edge count
        for fid in 1:n_faces, lid in 1:elem_dim
            he = half_edges[fid][lid]
            twin_idx = vert_he_map[he.dest, he.origin]
            if v2e[he.origin] == (0,0)
                v2e[he.origin] = (fid,lid)
            end
            if twin_idx != (0,0)
                e2e[fid,lid] = twin_idx
            else
                n_boundary += 1
                bhe = HalfEdge(he.dest, he.origin)
                if n_boundary <= size(half_edges,1)
                    push!(half_edges[n_boundary], bhe)
                else
                    push!(half_edges, [bhe])
                end
                boundary_idx = (n_boundary,0)
                e2e[fid,lid] = boundary_idx
                push!(b2e, (fid,lid))
                v2e[he.dest] = boundary_idx
            end
        end

        new(n_vertices, n_faces, n_boundary, vertices, half_edges, ec, v2e, e2e, b2e)
    end
end

# Uniform volume mesh structure
# TODO: Support non-uniform mesh
mutable struct VolumeMesh{T<:Real}
    ### Attributes
    # Number of vertices in the mesh
    n_vertices::Int64
    # Number of cells in the mesh
    n_cells::Int64
    # Array of vertices in the mesh
    vertices::Array{Vertex{3,T}}
    # Array of anchor half-faces in the mesh
    #   (indices are (cell ID, local ID, anchor ID)
    half_faces::Array{Array{Array{HalfFace}}}
    # Element connectivity array
    ec::Array{Int64, 2}
    # Map of each vertex to AHF anchored at the vertex
    #   & border vertex to boundary AHF
    v2f::Array{Tuple{Int64,Int64,Int64}}
    # Map of each internal half-face with anchor ID 1 to its twin AHF
    f2f::Array{Tuple{Int64,Int64,Int64},2}
    # Map of each boundary half-face with anchor ID 1 to its twin AHF
    b2f::Array{Tuple{Int64,Int64,Int64}}

    ### Constructor & check mesh sanity
    function VolumeMesh(
        vertices::Array{Vertex{3,T}},
        ec::Array{Int64, 2} )
    	# populate number of elements
        n_vertices = size(vertices)
        n_cells = size(cells)
        elem_dim = size(ec,2)

        # Hardcode half-face topology per cell
        cell_hf_topo = []
        if elem_dim == 4 # tetrahedron
            cell_hf_topo = [1 2 3; 1 4 2; 2 4 3; 1 3 4]
        elseif elem_dim == 6 # hexahedron
            cell_hf_topo = [1 2 3 4; 1 6 7 2; 2 7 8 3; 3 8 5 4; 1 4 5 6; 5 6 7 8]
        end
        face_size = size(cell_hf_topo,2)
        
        # Populate half-face array
        vert_hf_map = Dict(Array{Int64},HalfFace())
        half_faces = Vector{Vector{Vector{HalfFace}}}()
        for cid in 1:n_cells
            cell_half_faces = Vector{Vector{HalfFace}}()
            for lid in 1:elem_dim
                half_face_anchors = Array{HalfFace}(undef, face_size)
                hf_ids = [ec[cid,i] for i in cell_hf_topo[lid,:]]
                for aid in 1:face_size
                    hf_ids_shift = circshift(hf_ids, aid)
                    hf = HalfFace(hf_ids_shift, aid)
                    half_face_anchors[aid] = hf
                    vert_hf_map[hf_ids_shift] = hf
                end
                push!(cell_half_faces, half_face_anchors)
            end
            push!(half_faces, cell_half_faces)
        end

        # Populate half-face - twin half-face map
        f2f = Tuple{Int64,Int64,Int64}[(0,0,0) for i in 1:n_faces, j in 1:elem_dim]
        v2f = Tuple{Int64,Int64,Int64}[(0,0,0) for i in 1:n_vertices]
        b2e = Tuple{Int64,Int64,Int64}[]
        n_boundary = 0 # Boundary half-face count
        for cid in 1:n_cells, lid in 1:elem_dim
            hf = half_faces[cid][lid]
            twin_idx = circshift(get(vert_hf_map, reverse(hf.vertices), undef), 1)
            if v2f[hf.vertex[hf.anchor]] == (0,0,0)
            	v2f[hf.vertex[hf.anchor]] = (cid, lid, 1)
            end 
            if twin_idx != undef
            else
            end
        end
    end
    
end
