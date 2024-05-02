from calendar import c
from altair import layer
import gdsfactory as gf
from networkx import center

core_layer = (1, 0)


def phase_shifter(l, wwg, nbends, ubend):
    c = gf.Component()
    ls = l-2*ubend.size[0]
    s0 = c << gf.path.extrude(gf.path.straight(
        length=ls), width=wwg, layer=core_layer)
    s = s0

    lu = ubend.size[0]
    sin = c << gf.path.extrude(gf.path.straight(
        length=lu), width=wwg, layer=core_layer)
    sin.connect("o2", s.ports["o1"])

    for i in range(1, nbends):
        u_ = (c << ubend)
        p1 = "o1" if i % 2 == 1 else "o2"
        p2 = "o2" if i % 2 == 1 else "o1"
        u_.connect(p1, s.ports["o2"])
        s_ = (c << gf.path.extrude(gf.path.straight(
            length=ls), width=wwg, layer=core_layer))
        s_.connect("o1", u_.ports[p2])
        s = s_

        pn_layer = (71, 0) if i % 2 == 1 else (72, 0)
        pn = c << gf.components.rectangle(
            size=(l, ubend.size[1]-wwg), layer=pn_layer)
        pn.xmin = c.xmin
        pn.ymin = s.center[1]

    sout = c << gf.path.extrude(gf.path.straight(
        length=lu), width=wwg, layer=core_layer)
    sout.connect("o1", s.ports["o2"])

    heater_margin = 0.1
    center = c.center
    heater = c << gf.components.rectangle(
        size=(l+2*heater_margin, c.size[1]+2*heater_margin), layer=(61, 0))
    heater.center = center
    c.show()
    return c


def ubend_resonator(l, wwg, ubend):
    c = gf.Component()
    arm = gf.Component()

    r = 1
    lc = 8
    ls = (l-lc-4*r-2*ubend.size[0])/2

    P = gf.Path()
    P += gf.path.straight(length=ls)  # Straight section
    P += gf.path.euler(radius=r, angle=-45, use_eff=True)
    P += gf.path.euler(radius=r, angle=45, use_eff=True)
    P += gf.path.straight(length=lc/2)  # Straight section
    gf.components.bend_s((4, 1), 32, wwg, True)

    coupler_margin = .15
    quarter = gf.path.extrude(P, width=wwg, layer=core_layer)

    sw = arm << quarter
    sw.ymin = +wwg/2+coupler_margin
    sw.xmax = 0

    se = arm << quarter
    se = se.mirror()
    se.connect("o2", quarter.ports["o2"])

    arm.add_port("o1", port=sw.ports["o1"])
    arm.add_port("o2", port=se.ports["o1"])
    c << arm

    ubend1 = c << ubend
    ubend2 = c << ubend
    ubend1.connect("o1", arm.ports["o1"])
    ubend2.connect("o2", arm.ports["o2"])
    (c << arm).connect("o2", ubend1.ports["o2"])

    leg1 = c << gf.components.straight(length=l, width=wwg, layer=core_layer)
    # c.add_port(name="o1", center=(0, 0), width=wwg, orientation=0)
    leg1.center = (0, 0)
    y = c.ymax+coupler_margin
    leg2 = c << gf.components.straight(length=l, width=wwg, layer=core_layer)
    leg2.ymin = y
    c.show()
    return c
