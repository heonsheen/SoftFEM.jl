include("../src/Geometry.jl")

import Makie
import AbstractPlotting
import GeometryTypes
GT = GeometryTypes

#=
nx = 3
ny = 3
vertices = Array{Vertex}(undef, nx * ny)
for i in 0:ny-1, j in 0:nx-1
    vertices[i * nx + j + 1] = Vertex(Array{Float64}([2.0 - i, j, i]))
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
surf_mesh = extract_surface(volume_mesh)

scene = Makie.Scene()
node = Makie.Node(0.0)

vts = [v_i.x[j] for v_i in surf_mesh.vertices, j = 1:3]
s1 = Makie.mesh!(scene, vts, surf_mesh.ec, color = :blue, shading = false, show_axis = false)[end]
s2 = Makie.wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3, show_axis = false)[end]
Makie.display(scene)

#=
VT = vertextype(GT.GLNormalMesh)
FT = facetype(GT.GLNormalMesh)
=#

N = 50
#Makie.record(scene, "output/test.mp4", 1:N) do i
for i = 1:N
    volume_mesh.vertices[1].x[3] += 0.01
    volume_mesh.vertices[7].x[3] -= 0.01
    
    surf_mesh = extract_surface(volume_mesh)

    vts = [v_i.x[j] for v_i in surf_mesh.vertices, j = 1:3]
    #=
    vs = [VT(v.x[1], v.x[2], v.x[3]) for v in surf_mesh.vertices]
    fs = [FT(surf_mesh.ec[fi,1], surf_mesh.ec[fi,2], surf_mesh.ec[fi,3]) for fi in 1:size(surf_mesh.ec,1)]
    glmesh = GT.GLNormalMesh(vs, fs)
    s[1] = glmesh
    =#
    s1[1] = vts
    #AbstractPlotting.force_update!()
    sleep(1/24)
end
