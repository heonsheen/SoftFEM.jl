include("src/Geometry.jl")
include("src/CGTriObject.jl")
include("src/CGTetObject.jl")
include("src/LinearElasticMaterial.jl")
include("src/NeohookeanMaterial.jl")
include("src/BackwardEuler.jl")

import Makie
#import AbstractPlotting
import GeometryTypes
using UnicodePlots
#using IterativeSolvers

GT = GeometryTypes

#=
nx = 3
ny = 3
vertices = Array{Vertex}(undef, nx * ny)
for i in 0:ny-1, j in 0:nx-1
    vertices[i * nx + j + 1] = Vertex(Array{Float64}([2.0 - i, j]))
end

ec = [1 5 4;
      1 2 5;
      2 3 5;
      3 6 5;
      6 9 5;
      8 5 9;
      7 5 8;
      5 7 4]

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

ec = [2 3 7 6;
      2 6 7 5;
      3 4 7 6;
      4 5 7 6;
      2 6 1 3;
      5 6 1 2;
      6 4 1 3;
      6 5 1 4]

mesh = VolumeMesh(vertices, ec)
surf_mesh = extract_surface(mesh)

mp = Dict{String,Float64}(
    "E" => 1.0,
    "nu" => 0.35
)
mat = NeohookeanMaterial(mp, 0.01)

obj = CGTetObject(mesh, mat)

n_steps = 500
N = obj.N
dim = obj.dim

fixed = zeros(Bool, N*dim)
for i in [19, 20, 21]
      fixed[i] = true
end

dt = 0.01
g = repeat([0.0; 0.0; -9.81], N)
#u = zeros(N*dim*2)

#limits = Makie.IRect(-5, -5, 10, 10)
scene = Makie.Scene()
node = Makie.Node(0.0)

vts = [v_i.x[j] for v_i in surf_mesh.vertices, j = 1:3]
s1 = Makie.mesh!(scene, vts, surf_mesh.ec, color = :blue, shading = false, show_axis = false)[end]
s2 = Makie.wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3, show_axis = false)[end]
Makie.display(scene)

for timestep = 1:n_steps
#while true
      u = [obj.x - obj.X; obj.v]
      u_new = backward_euler(u, obj, dt, fixed, g)
      dx = u_new[1:N*dim]
      v = u_new[N*dim+1:end]

      update_mesh(mesh, obj)
      surf_mesh = extract_surface(mesh)
      vts = [v_i.x[j] for v_i in surf_mesh.vertices, j = 1:3]
      s1[1] = vts

      #u = u_new
      sleep(1/24)
end
