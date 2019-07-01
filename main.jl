include("src/Geometry.jl")
include("src/CGTriObject.jl")
include("src/LinearElasticMaterial.jl")
include("src/BackwardEuler.jl")

import Makie
#import AbstractPlotting
import GeometryTypes
using UnicodePlots
#using IterativeSolvers

GT = GeometryTypes


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
#=
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
surf_mesh = extract_surface(volume_mesh)
=#

mp = Dict{String,Float64}(
    "E" => 1.0,
    "nu" => 0.35
)
mat = LinearElasticMaterial(mp, 0.01)

obj = CGTriObject(mesh, mat)

n_steps = 100
N = obj.N
dim = obj.dim

fixed = zeros(Bool, N*dim)
for i in [5,6,11,12,17,18]
      fixed[i] = true
end

dt = 0.01
g = repeat([0.0; -9.81], N)
#u = zeros(N*dim*2)

limits = Makie.IRect(-5, -5, 10, 10)
scene = Makie.Scene(limits = limits)
node = Makie.Node(0.0)

vts = [v_i.x[j] for v_i in mesh.vertices, j = 1:2]
s1 = Makie.mesh!(scene, vts, mesh.ec, color = :blue, shading = false, show_axis = false)[end]
s2 = Makie.wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3, show_axis = false)[end]
Makie.display(scene)

for timestep = 1:n_steps
      u = [obj.x - obj.X; obj.v]
      u_new = backward_euler(u, obj, dt, fixed, g)
      dx = u_new[1:N*dim]
      v = u_new[N*dim+1:end]

      update_mesh(mesh, obj)
      vts = [v_i.x[j] for v_i in mesh.vertices, j = 1:2]
      s1[1] = vts

      #u = u_new
      sleep(1/24)
end
