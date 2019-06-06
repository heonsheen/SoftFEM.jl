# Simple definitions for basic geometry structure

using StaticArrays
using LinearAlgebra

# N-dimensional Vertex data with abstract Real type T
mutable struct Vertex{N<:Int64, T<:Real}
### Attributes
    # Positions of the vertex
    x::SVector{N,T}

### Constructors
    function Vertex{N,T}(x::SVector{N,T}) where {N<:Int64, T<:Real}
    	if N < 2
    		error("Dimension cannot be less than 2")
    	elseif N > 3
    		error("Dimension cannot be larger than 3")
    	end
        new(x)
    end
end

# Half-edge connectivity data
mutable struct HalfEdge
### Attributes	
	# Index of the half-edge after
	next::Int64
	# Index of the half-edge before
	prev::Int64
	# Index of the twin half-edge
	twin::Int64
	# Index of the Face that the half-edge belongs to
	# (if the face is a hole, index -1)
	face::Int64
	# Index of the vertex the half-edge originates from
	origin::Int64
	# NOTE: Destination vertex index givin by 
	# 		dest = he.next.origin

### Constructors
	function HalfEdge(
		next::Int64,
		prev::Int64,
		twin::Int64,
		face::Int64,
		origin::Int64)
		if next < 0 || pref < 0 || twin < 0 || origin < 0
			error("Invalid indexing.")
		end
		new(next, prev, twin, face, origin)
	end
end	

# Face connectivity data for triangle topology
mutable struct Face2D
### Attributes
	# The index of a half-edge belonging to this face
	half_edge::Int64
	#=
    # Vertices Topology
    topo::Tuple{Int64,Int64,Int64}
    # Half-edge Topology
    halfEdges::Tuple{Int64,Int64,Int64}
	=#

### Constructor
	function Face2D(half_edge::Int64)
		if half_edge < 0
			error("invalid indexing")
		end
		new(half_edge)
	end
	#=
    function Face(
    	topo::Tuple{Int64,Int64,Int64},
        halfEdges::Tuple{Node2d{T},Node2d{T},Node2d{T}} )
    	if topo[1] <= 0 || topo[2] <= 0 || topo[3] <= 0
    		error("Invalid vertex topology indexing")
    	elseif halfEdges[1] <= 0 || halfEdges[2] <= 0 || halfEdges[3] <= 0
    		error("Invalid half-edge topology indexing")
    	end
        new(topo, half-edges)
    end
    =#
end

# Half-Face connectivity data for volumetric mesh topology
mutable struct HalfFace
### Attributes
	# The index of a half-edge belonging to this face
	half_edge::Int64
	# The index of the cell the half-face belongs to
	cell::Int64

### Constructor
	function HalfFace(half_edge::Int64, cell::Int64)
		if half_edge < 0 || cell < 0
			error("invalid indexing")
		end
		new(half_edge, cell)
	end
end

# Cell connectivity data 
mutable struct Cell
### Attributes
	# The index of a half-face beloning to this cell
	half_face::Int64

### Constructor
	function Cell(half_face::Int64)
		if half_face < 0
			error("invalid indexing")
		end
		new(half_face)
	end
end

# Triangle mesh structure with dimensionality N
mutable struct Trimesh{N<:Int64, T<:Real}
### Attributes
	# Number of vertices in the mesh
    n_verts::Int64
    # Number of half-edges in the mesh
    n_half_edges::Int64
    # Number of faces in the mesh
    n_faces::Int64
    # Array of vertices in the mesh
    vertices::Array{Vertex{N,T}}
    # Array of half-edges in the mesh
    half_edges::Array{HalfEdge}
    # Array of faces in the mesh
    faces::Array{Face}
    # Array of internal node indices
    internal_nodes_topo::Array{Int64}
    # Array of boundary node indices
    boundary_nodes_topo::Array{Int64}
    # Adjacency list for vertices
    adjacency_list::Array{Int64,2}

### Constructor & check mesh sanity
    function Trimesh(
        vertices::Array{Vertex{N,T}},
        half_edges::Array{HalfEdge},
        faces::Array{Face} ) where {N<:Int64, T<:Real}
    	# populate number of elements
        n_verts = size(vertices)
        n_half_edges = size(half_edges)
        n_faces = size(faces)

        ### Check mesh sanity
        # Check face.half_edge.face == face
       	for fid in 1:n_faces
       		face = faces[fid]
       		he = half_edges[face.half_edge]
       		if he.face != fid
       			error("face " * fid * ": face.half_edge.face != face")
       		end
   		end

   		# Check if half-edges in a face construct a loop
   		for fid in 1:n_faces
   			
   		end
    end
    
end
