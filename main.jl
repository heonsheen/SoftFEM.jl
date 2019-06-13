include("Geometry.jl")
#=
nx = 3
ny = 3
vertices = Array{Vertex}(undef, nx * ny)
for i in 0:ny-1, j in 0:nx-1
    vertices[i * nx + j + 1] = Vertex(Array{Float64}([2.0 - i, j]))
end

ec = [1 4 5;
      1 5 2;
      2 5 3;
      3 5 6;
      6 5 9;
      8 9 5;
      7 8 5;
      5 4 7]

mesh = Mesh(vertices, ec)
=#

points = [0 0 -1;
		  sqrt(2) -sqrt(2) 0;
		  sqrt(2) sqrt(2) 0;
		  -sqrt(2) sqrt(2) 0;
		  -sqrt(2) -sqrt(2) 0;
		  0 0 0;
		  0 0 1]

n_points = size(points, 1)
vertices = Array{Vertex}(undef, n_points)
for i in 1:n_points
	vertices[i] = Vertex(Array{Float64}(points[i,:]))
end

ec = [2 3 6 7;
      2 6 5 7;
      3 4 6 7;
      4 5 6 7;
      2 6 3 1;
      5 6 2 1;
      6 4 3 1;
      6 5 4 1]

 volume_mesh = VolumeMesh(vertices, ec)