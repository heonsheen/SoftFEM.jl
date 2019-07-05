include("src/Geometry.jl")
include("src/CGTriObject.jl")
include("src/DGTriObject.jl")
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


nx = 3
ny = 3
vertices = Array{Vertex}(undef, nx * ny)
for i in 0:ny-1, j in 0:nx-1
    vertices[i * nx + j + 1] = Vertex(Array{Float64}([-1.0+j, -1.0+i]))
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

fixed = zeros(Bool, nx*ny*2)
for i in [13, 14, 15, 16, 17, 18]
      fixed[i] = true
end

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

ec = [2 3 7 6;
      2 6 7 5;
      3 4 7 6;
      4 5 7 6;
      2 6 1 3;
      5 6 1 2;
      6 4 1 3;
      6 5 1 4]
=#
#=
nx = 10
ny = 10
nz = 10
dx = 1.0 / (nx-1)
dy = 1.0 / (ny-1)
dz = 1.0 / (nz-1)

fixed = zeros(Bool, nx*ny*nz*3)

points = zeros(Float64, (nx)*(ny)*(nz), 3)
for k = 0:nz-1, j = 0:ny-1, i = 0:nx-1
      x = i + j*nx + k*nx*ny + 1
      points[x, :] = [-0.5 + i*dx, -0.5 + j*dy, -0.5 + k*dz]
      if k == nz-1
            for d = 1:3
                  fixed[3*(x-1)+d] = true
            end
      end
end
n_points = size(points, 1)
vertices = Array{Vertex}(undef, n_points)
for i in 1:n_points
	vertices[i] = Vertex(Array{Float64}(points[i,:]))
end

ec = zeros(Int64, 5*(nx-1)*(ny-1)*(nz-1), 4)
for k = 0:nz-2, j = 0:ny-2, i = 0:nx-2
      cell_id = i + j*(nx-1) + k*(nx-1)*(ny-1) + 1
      x = i + j*nx + k*nx*ny + 1
      if mod(cell_id,2) == 1
            ec[5*(cell_id-1)+1,:] = [x, x+1, x+nx*ny+1, x+nx+1]
            ec[5*(cell_id-1)+2,:] = [x, x+nx*ny+1, x+nx*(ny+1), x+nx+1]
            ec[5*(cell_id-1)+3,:] = [x, x+nx*ny+1, x+ny*ny, x+nx*(ny+1)]
            ec[5*(cell_id-1)+4,:] = [x, x+nx*(ny+1), x+nx, x+nx+1]
            ec[5*(cell_id-1)+5,:] = [x+nx+1, x+nx*(ny+1), x+nx*(ny+1)+1, x+nx*ny+1]
      else
            ec[5*(cell_id-1)+1,:] = [x, x+1, x+nx*ny, x+nx]
            ec[5*(cell_id-1)+2,:] = [x+1, x+nx+1, x+nx*(ny+1)+1, x+nx]
            ec[5*(cell_id-1)+3,:] = [x+1, x+nx*ny, x+nx, x+nx*(ny+1)+1]
            ec[5*(cell_id-1)+4,:] = [x+1, x+nx*(ny+1)+1, x+nx*ny+1, x+nx*ny]
            ec[5*(cell_id-1)+5,:] = [x+nx, x+nx*ny, x+nx*(ny+1), x+nx*(ny+1)+1]
      end
end
=#
#=
mesh = VolumeMesh(vertices, ec)
surf_mesh = extract_surface(mesh)
=#
mp = Dict{String,Float64}(
    "E" => 0.5,
    "nu" => 0.35
)
mat = NeohookeanMaterial(mp, [0.1, 0.1], 0.01)

obj = DGTriObject(mesh, mat)

n_steps = 100
N = obj.N
dim = obj.dim

dt = 0.01
#g = repeat([0.0; 0.0; -9.81], N)
g = repeat([0.0; -9.81], N)
#u = zeros(N*dim*2)

#limits = Makie.IRect(-5, -5, 10, 10)
scene = Makie.Scene()
node = Makie.Node(0.0)
#=
vts = [v_i.x[j] for v_i in surf_mesh.vertices, j = 1:3]
s1 = Makie.mesh!(scene, vts, surf_mesh.ec, color = :blue, shading = false, show_axis = false)[end]
=#
vts = [v_i.x[j] for v_i in mesh.vertices, j = 1:2]
s1 = Makie.mesh!(scene, vts, mesh.ec, color = :blue, shading = false, show_axis = false)[end]

s2 = Makie.wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3, show_axis = false)[end]
Makie.display(scene)

#for timestep = 1:n_steps
#while true
Makie.record(scene, "results/video.mp4", 1:n_steps) do timestep
      u = [obj.x - obj.X; obj.v]
      u_new = backward_euler(u, obj, dt, fixed, g)
      dx = u_new[1:N*dim]
      v = u_new[N*dim+1:end]

      update_mesh(mesh, obj)
      #surf_mesh = extract_surface(mesh)
      #vts = [v_i.x[j] for v_i in surf_mesh.vertices, j = 1:3]
      vts = [v_i.x[j] for v_i in mesh.vertices, j = 1:2]
      s1[1] = vts

      #u = u_new
      println("timestep ", timestep)
      sleep(1/120)
end
