# from luminescent import finetune, apply_design, load_component, load_solution
# # sol = finetune(2)
# sol = load_solution(study="inverse_design")
# c = load_component(study="inverse_design")
# c = apply_design(c, sol)
import gdsfactory as gf

c = gf.Component()
e0 = c << gf.components.ellipse(layer=(1, 0))
e1 = c << gf.components.ellipse(layer=(1, 0))
e2 = c << gf.components.ellipse(layer=(1, 0))
e3 = c << gf.components.ellipse(layer=(1, 0))
e4 = c << gf.components.ellipse(layer=(1, 0))
e5 = c << gf.components.ellipse(layer=(1, 0))

e1.drotate(15 * 1)
e2.drotate(15 * 2)
e3.drotate(15 * 3)
e4.drotate(15 * 4)
e5.drotate(15 * 5)

c
polygons = c.get_polygons(merge=True)
p = polygons[1][0]
r = gf.Component()
r.add_polygon(p, layer=(1, 0))
r.show()
