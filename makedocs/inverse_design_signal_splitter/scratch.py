import json
import os

dir = "_temp"
p = os.path.join(dir, "sol.json")
# sol = bson.loads(p, "rb").read())["sol"]
sol = json.loads(open(p).read())
sol["sparams"] = {k: (v[0]+1j*v[1])
                  for k, v in sol["sparams"].items()}
# sol["designs"] = [np.array(d) for d in sol["designs"]]
raise ValueError("stop here")
scene = c.to_3d(layer_stack=LAYER_STACK)
# scene.show()

g = sum(scene.geometry.values())
s = g.section(plane_origin=g.centroid,
              plane_normal=[0,  1, 0])
p, _ = s.to_planar()
p = clip_by_rect(GeometryCollection(p.polygons_full), -5, -5, 5, 5)

fig = plt.figure(1, )

# 3: invalid polygon, ring touch along a line
ax = fig.add_subplot(121)
plot_polygon(p, ax=ax)
# plt.show()

# img = rasterio.features.rasterize([p], out_shape=)
plt.imshow(img)

pdk = get_generic_pdk()
pdk.activate()

component = gf.components.coupler_ring()
component.plot()

# define a mapping of pdk material names to tidy3d medium objects
mapping = {
    "si": td.Medium(name="Si", permittivity=3.47**2),
    "sio2": td.Medium(name="SiO2", permittivity=1.47**2),
}

# setup the tidy3d component
c = gt.Tidy3DComponent(
    component=component,
    layer_stack=LAYER_STACK,
    material_mapping=mapping,
    pad_xy_inner=2.0,
    pad_xy_outer=2.0,
    pad_z_inner=0,
    pad_z_outer=0,
    extend_ports=2.0,
)
c.slice_stack((1))

# plot the component and the layerstack
fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(ncols=2, nrows=3, width_ratios=(3, 1))
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[1, 0])
ax2 = fig.add_subplot(gs[2, 0])
axl = fig.add_subplot(gs[1, 1])
c.plot_slice(x="core", ax=ax0)
c.plot_slice(y="core", ax=ax1)
c.plot_slice(z="core", ax=ax2)
axl.legend(*ax0.get_legend_handles_labels(), loc="center")
axl.axis("off")
plt.show()
