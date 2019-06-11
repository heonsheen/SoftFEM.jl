# Simple definitions for basic geometry structure

using StaticArrays
using LinearAlgebra

# N-dimensional Vertex data with abstract Real type T
mutable struct Vertex{N<:Int64, T<:Real}
### Attributes
    # Positions of the vertex
    x::SVector{N,T}
    # Index of a half-edge going out of this vertex
    half_edge::Int64

### Constructors
    function Vertex{N,T}(x::SVector{N,T}, half_edge::Int64) where {N<:Int64, T<:Real}
    	@assert N >= 2 "Dimension cannot be less than 2"
    	@assert N <= 3 "Dimension cannot be larger than 3"
    	
        new(x, half_edge)
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
mutable struct Face
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
	function Face(half_edge::Int64)
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
	# Index of the half-edge after
	next::Int64
	# Index of the half-edge before
	prev::Int64
	# Index of the twin half-edge
	twin::Int64
	# The index of a half-edge belonging to this face
	half_edge::Int64
	# The index of the cell the half-face belongs to
	# note: if the cell is a hole, then cell = -1
	cell::Int64

### Constructor
	function HalfFace(
		next::Int64,
		prev::Int64,
		twin::Int64,
		half_edge::Int64, 
		cell::Int64)
		if next < 0 || prev < 0 || twin < 0 || half_edge < 0 || cell < 0
			error("invalid indexing")
		end
		new(next, prev, twin, half_edge, cell)
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
    n_vertices::Int64
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
    # Array of boundary half-edge indices
    boundary_edges_topo::Array{Int64}
    # Adjacency list for vertices
    adjacency_list::Array{Int64,2}

### Constructor & check mesh sanity
    function Trimesh(
        vertices::Array{Vertex{N,T}},
        half_edges::Array{HalfEdge},
        faces::Array{Face} ) where {N<:Int64, T<:Real}
    	# populate number of elements
        n_vertices = size(vertices)
        n_half_edges = size(half_edges)
        n_faces = size(faces)

        ### Check mesh sanity
        # Check face.half_edge.face == face
       	for fid in 1:n_faces
       		face = faces[fid]
       		he = half_edges[face.half_edge]
       		@assert he.face == fid ("face " * fid * ": face.half_edge.face != face")
   		end

   		# Check if half-edges in a face construct a loop
   		max_face_half_edges = 100
   		for fid in 1:n_faces
   			num_face_half_edges = 0
   			face = faces[fid]
   			he = face.half_edge
   			he_end = face.half_edge
   			while true
   				num_face_half_edges += 1
   				@assert num_face_half_edges < max_face_half_edges ("face " * fid * " has too many half-edges")
   				he = half_edges[he].next
   				he == he_end && break
   			end
   		end

   		# Check mutual relationship between half-edges
   		for heid in 1:n_half_edges
   			he = half_edges[heid]
   			@assert half_edges[he.twin].twin == heid ("half-edge " * heid * " != he.twin.twin")
   			@assert half_edges[he.next].prev == heid ("half-edge " * heid * " != he.next.prev")
   			@assert half_edges[he.prev].next == heid ("half-edge " * heid * " != he.prev.next")
   			@assert he.twin != heid ("half-edge " * heid * " == he.twin")
   			@assert he.prev != heid ("half-edge " * heid * " == he.prev")
   			@assert he.next != heid ("half-edge " * heid * " == he.next")
   		end

   		# Vertices sanity check
   		max_outgoing_half_edges = 100
   		for vid in 1:n_vertices
   			n_traversed_edges = 0
   			vert = vertices[vid]
   			is_vertex_visited = Zeros(Bool, n_vertices)
   			boundary_flag = false

   			# Vertex ring iteration
   			ring_iter = half_edges[vert.half_edge].twin
   			ri_0 = ring_iter
   			while true
   				other_vertex = half_edges[ring_iter].origin
   				@assert half_edges[ring_iter].dest == vid ("ring_iter.dest != dest")
   				@assert num_traversed_edges < max_outgoing_half_edges ("Too many half-edges in ring")
   				@assert is_vertex_visited[other_vertex] == false ("More than one edge between two vertices")
   				is_vertex_visited[other_vertex] = true
   				num_traversed_edges += 1

   				# populate adjacency list
   				adjacency_list[vid, other_vertex] = 1

   				ring_iter = half_edges[half_edges[ring_iter].next].twin
   				ring_iter == ri_0 && break
   			end
   		end

   		# Identify boundary half-edges
   		for heid in 1:n_half_edges
   			if half_edges[heid].face == -1
   				bhe = heid
   				do
   					push!(boundary_edges_topo, bhe)
   					bhe = half_edges[bhe].next
	   				bhe == heid && break
	   			end
   				break
   			end
   		end
    end
end


# Tet mesh structure
mutable struct Tetmesh{T<:Real}
### Attributes
	# Number of vertices in the mesh
    n_vertices::Int64
    # Number of half-edges in the mesh
    n_half_edges::Int64
    # Number of faces in the mesh
    n_half_faces::Int64
    # Number of cells in the mesh
    n_cells::Int64
    # Array of vertices in the mesh
    vertices::Array{Vertex{3,T}}
    # Array of half-edges in the mesh
    half_edges::Array{HalfEdge}
    # Array of faces in the mesh
    half_faces::Array{HalfFace}
    # Array of cells in the mesh
    cells::Array{Cell}
    # Array of internal node indices
    internal_nodes_topo::Array{Int64}
    # Array of boundary node indices
    boundary_nodes_topo::Array{Int64}
    # Adjacency list for vertices
    adjacency_list::Array{Int64,2}

### Constructor & check mesh sanity
    function Tetmesh(
        vertices::Array{Vertex{3,T}},
        half_edges::Array{HalfEdge},
        half_faces::Array{HalfFace},
        cells::Array{Cell} ) where {T<:Real}
    	# populate number of elements
        n_vertices = size(vertices)
        n_half_edges = size(half_edges)
        n_half_faces = size(half_faces)
        n_cells = size(cells)

        ### Check mesh sanity
        # Check half_face.half_edge.face == half_face
       	for hfid in 1:n_half_faces
       		half_face = half_faces[hfid]
       		he = half_edges[half_face.half_edge]
       		@assert he.face == hfid ("face " * fid * ": face.half_edge.face != face")
   		end

   		# Check if half-edges in a half-face construct a loop
   		max_half_face_half_edges = 100
   		for hfid in 1:n_half_faces
   			num_half_face_half_edges = 0
   			half_face = half_faces[hfid]
   			he = half_face.half_edge
   			he_end = half_face.half_edge
   			while true
   				num_half_face_half_edges += 1
   				@assert num_half_face_half_edges < max_half_face_half_edges ("half_face " * hfid * " has too many half-edges")
   				he = half_edges[he].next
   				he == he_end && break
   			end
   		end

   		# Check if half-faces in a cell construct a loop
   		max_cell_half_faces = 100
   		for cid in 1:n_cells
   			num_cell_half_faces = 0
   			cell = cells[cid]
   			hf = cell.half_face
   			hf_end = cell.half_face
   			while true
   				num_cell_half_faces += 1
   				@assert num_cell_half_faces < max_cell_half_faces ("cell " * cid * " has too many half-faces")
   				hf = half_faces[hf].next
   				hf == hf_end && break
   			end
   		end

   		# Check mutual relationship between half-edges
   		for heid in 1:n_half_edges
   			he = half_edges[heid]
   			@assert half_edges[he.twin].twin == heid ("half-edge " * heid * " != he.twin.twin")
   			@assert half_edges[he.next].prev == heid ("half-edge " * heid * " != he.next.prev")
   			@assert half_edges[he.prev].next == heid ("half-edge " * heid * " != he.prev.next")
   			@assert he.twin != heid ("half-edge " * heid * " == he.twin")
   			@assert he.prev != heid ("half-edge " * heid * " == he.prev")
   			@assert he.next != heid ("half-edge " * heid * " == he.next")
   		end

   		# Check mutual relationship between half-faces
   		for hfid in 1:n_half_faces
   			hf = half_faces[hfid]
   			@assert half_faces[hf.twin].twin == hfid ("half-face " * hfid * " != hf.twin.twin")
   			@assert half_faces[hf.next].prev == hfid ("half-face " * hfid * " != hf.next.prev")
   			@assert half_faces[hf.prev].next == hfid ("half-face " * hfid * " != hf.prev.next")
   			@assert hf.twin != heid ("half-face " * hfid * " == hf.twin")
   			@assert hf.prev != heid ("half-face " * hfid * " == hf.prev")
   			@assert hf.next != heid ("half-face " * hfid * " == hf.next")
   		end

   		# Vertices sanity check
   		max_outgoing_half_edges = 100
   		for vid in 1:n_vertices
   			n_traversed_edges = 0
   			vert = vertices[vid]
   			is_vertex_visited = Zeros(Bool, n_vertices)
   			boundary_flag = false

   			# Vertex ring iteration
   			ring_iter = half_edges[vert.half_edge].twin
   			ri_0 = ring_iter
   			while true
   				other_vertex = half_edges[ring_iter].origin
   				@assert half_edges[ring_iter].dest == vid ("ring_iter.dest != dest")
   				@assert num_traversed_edges < max_outgoing_half_edges ("Too many half-edges in ring")
   				@assert is_vertex_visited[other_vertex] == false ("More than one edge between two vertices")
   				is_vertex_visited[other_vertex] = true
   				num_traversed_edges += 1

   				# Identify boundary nodes
   				ring_iter_twin = half_edges[ring_iter].twin
   				if half_edges[ring_iter].face == -1 || half_edges[ring_iter_twin].face == -1
   					boundary_flag = true
   				end

   				# populate adjacency list
   				adjacency_list[vid, other_vertex] = 1

   				ring_iter = half_edges[half_edges[ring_iter].next].twin
   				ring_iter == ri_0 && break
   			end

   			if boundary_flag
   				push!(boundary_nodes_topo, vid)
   			else
   				push!(internal_nodes_topo, vid)
   			end
   		end
    end
    
end

