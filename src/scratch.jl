using Meshes, GeoIO, GLMakie
using Meshes: Point, boundingbox

mesh = getfield(GeoIO.load("slab.obj"), :domain)
mesh_bbox = boundingbox(mesh)

function test(m)
    vol = map(Base.product(range.((0, -0.5, -0.5), (0.5, 0.5, 0.7), step=0.1)...)) do v
        Point(v) ∈ m && return 1
        0
    end
    GLMakie.volume(vol) |> display
end

test(mesh)
# test(mesh_bbox)