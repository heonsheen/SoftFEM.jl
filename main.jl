include("src/Geometry.jl")
include("src/CGTriObject.jl")
include("src/DGTriObject.jl")
include("src/CGTetObject.jl")
include("src/LinearElasticMaterial.jl")
include("src/NeohookeanMaterial.jl")
include("src/BackwardEuler.jl")
include("src/ForwardEuler.jl")
include("src/DGMixedIntegrator.jl")

import Makie
import AbstractPlotting
import GeometryTypes
using UnicodePlots
#using IterativeSolvers

GT = GeometryTypes

#=
nx = 3
ny = 3
vertices = Array{Vertex}(undef, nx * ny)
for j in 0:ny-1, i in 0:nx-1
    vertices[j * nx + i + 1] = Vertex(Array{Float64}([-1.0+i, -1.0+j]))
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
=#
#=
vertices = Array{Vertex}(undef, 4)
vertices[1] = Vertex(Array{Float64}([0.0, 1.0]))
vertices[2] = Vertex(Array{Float64}([-0.5, 0.0]))
vertices[3] = Vertex(Array{Float64}([0.5, 0.0]))
vertices[4] = Vertex(Array{Float64}([0.0, -1.0]))

ec = [1 2 3; 2 4 3]

mesh = Mesh(vertices, ec)
fixed = zeros(Bool, 8)
fixed[1] = true
fixed[2] = true
=#
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

nx = 5
ny = 5
dx = 1.0 / (nx-1)
dy = 1.0 / (ny-1)

fixed = zeros(Bool, nx*ny*2)

points = zeros(Float64, nx*ny, 2)
for j = 0:ny-1, i = 0:nx-1
      x = i + j*nx + 1
      points[x, :] = [-0.5 + i*dx, -0.5 + j*dy]
      if j == ny-1
            for d = 1:2
                  fixed[2*(x-1)+d] = true
            end
      end
end

n_points = size(points,1)
vertices = Array{Vertex}(undef, n_points)
for i in 1:n_points
      vertices[i] = Vertex(Array{Float64}(points[i,:]))
end

ec = zeros(Int64, 2*(nx-1)*(ny-1), 3)
for j = 0:ny-2, i = 0:nx-2
      cid = i + j * (nx-1) + 1
      x = i + j*nx + 1
      if mod(i+j,2) == 1
            ec[2*(cid-1)+1,:] = [x, x+1, x+nx]
            ec[2*(cid-1)+2,:] = [x+1, x+nx+1, x+nx]
      else
            ec[2*(cid-1)+1,:] = [x, x+1, x+nx+1]
            ec[2*(cid-1)+2,:] = [x, x+nx+1, x+nx]
      end
end

mesh = Mesh(vertices, ec)
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

mesh = VolumeMesh(vertices, ec)
surf_mesh = extract_surface(mesh)
=#
mp = Dict{String,Float64}(
    "E" => 1.0,
    "nu" => 0.35
)
mat = NeohookeanMaterial(mp, [0.0, 0.1], 0.05, 0.5, true)

obj = DGTriObject(mesh, mat)
dg_mesh = get_DG_mesh(obj)
dg_fixed = map_to_DG(obj, fixed)

n_steps = 200
N = obj.N
dim = obj.dim

dt = 0.01
#g = repeat([0.0; 0.0; -9.81], N)
g = repeat([0.0; -9.81], N)
#u = zeros(N*dim*2)
dg_g = map_to_DG(obj, g)

#limits = Makie.IRect(-5, -5, 10, 10)
scene = Makie.Scene(resolution = (750, 750))
node = Makie.Node(0.0)
#=
vts = [v_i.x[j] for v_i in surf_mesh.vertices, j = 1:3]
s1 = Makie.mesh!(scene, vts, surf_mesh.ec, color = :blue, shading = false, show_axis = false)[end]
=#
vts = [v_i.x[j] for v_i in dg_mesh.vertices, j = 1:2]
s1 = Makie.mesh!(scene, vts, dg_mesh.ec, color = :blue, 
                  shading = false, show_axis = false)[end]
s2 = Makie.wireframe!(scene[end][1], color = (:black, 0.6), 
                        linewidth = 3, show_axis = false)[end]

camera = Makie.cameracontrols(scene)
#=
Makie.@extractvalue camera (fov, near, projectiontype, lookat, eyeposition, upvector)
dir_vector = eyeposition - lookat
new_eyeposition = lookat + dir_vector * (1.0f0)
Makie.update_cam!(scene, new_eyeposition, lookat)
=#
Makie.update_cam!(scene, camera, 
                  GT.HyperRectangle{2,Float32}(Float32[-1.25, -1.75], Float32[2.5, 2.5]))
scene.center = false

#for timestep = 1:n_steps
#while true
Makie.record(scene, "results/video.mp4", 1:n_steps) do timestep
      u = [obj.x - obj.X; obj.v]
      u_new = dg_mixed_integrator(u, obj, dt, dg_fixed, dg_g)
      #u_new = backward_euler(u, obj, dt, dg_fixed, dg_g)
      dx = u_new[1:N*dim]
      v = u_new[N*dim+1:end]

      #update_mesh(mesh, obj)
      dg_mesh = get_DG_mesh(obj)
      #surf_mesh = extract_surface(mesh)
      #vts = [v_i.x[j] for v_i in surf_mesh.vertices, j = 1:3]
      vts = [v_i.x[j] for v_i in dg_mesh.vertices, j = 1:2]
      s1[1] = vts

      #u = u_new
      println("timestep ", timestep)
      sleep(1/120)
end
