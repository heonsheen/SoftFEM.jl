include("Geometry.jl")

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
