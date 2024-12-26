import gdsfactory as gf
from layers import *
core_layer = (1, 0)


def phase_shifter(l, wwg, nbends, ubend):
    c = gf.Component()
    ls = l-2*ubend.size[0]
    s0 = c << gf.path.extrude(gf.path.straight(
        length=ls), width=wwg, layer=core_layer)
    s = s0

    lu = ubend.size[0]
    wdope = ubend.size[1]-wwg
    sin = c << gf.path.extrude(gf.path.straight(
        length=lu), width=wwg, layer=core_layer)
    sin.connect("o2", s.ports["o1"])

    for i in range(nbends+2):
        if 0 < i < nbends+1:
            u_ = (c << ubend)
            p1 = "o1" if i % 2 == 1 else "o2"
            p2 = "o2" if i % 2 == 1 else "o1"
            u_.connect(p1, s.ports["o2"])
            s_ = (c << gf.path.extrude(gf.path.straight(
                length=ls), width=wwg, layer=core_layer))
            s_.connect("o1", u_.ports[p2])
            s = s_

        pn = "p" if i % 2 == 0 else "n"
        pn_layer = p_layer if pn == "p" else n_layer
        label = "p doping" if pn == "p" else "n doping"

        pn = c << gf.components.rectangle(size=(l, wdope), layer=pn_layer)
        pn.xmin = c.xmin
        if i < nbends+1:
            pn.ymin = s.center[1]
        else:
            pn.ymax = s.center[1]
        c.add_label(label, position=pn.center, layer=(501, 0))

    sout = c << gf.path.extrude(gf.path.straight(
        length=lu), width=wwg, layer=core_layer)
    sout.connect("o1", s.ports["o2"])

    heater_margin = 0.1
    center = c.center
    heater = c << gf.components.rectangle(
        size=(l+2*heater_margin, c.size[1]+2*heater_margin), layer=(61, 0))
    heater.center = center
    c.add_label("heater pad", position=heater.bbox_np()[1], layer=(502, 0))

    c.add_label("zigzag phase shifter", position=(
        c.center[0], c.ymax), layer=(500, 0))
    c.show()
    return c


def ubend_resonator(l, wwg, ubend):
    c = gf.Component()
    arm = gf.Component()
    quarter = gf.Component()

    r = 1
    common_left_pad_amount = 8
    lsb = 6
    ls = (l-common_left_pad_amount-2*lsb-2*ubend.size[0])/2

    su = quarter << gf.components.straight(
        length=ls, width=wwg, layer=core_layer)  # Straight section
    sc = quarter << gf.components.straight(
        length=common_left_pad_amount/2, width=wwg, layer=core_layer)  # Straight section
    sb = quarter << gf.components.bend_s(
        (lsb, 1), 32, gf.cross_section.strip(wwg), )
    sb.mirror(p1=(0, 0), p2=(1, 0))
    sb.connect("o2", sc.ports["o1"])
    su.connect("o2", sb.ports["o1"])
    quarter.add_port("o1", port=su.ports["o1"])
    quarter.add_port("o2", port=sc.ports["o2"])

    coupler_margin = .15

    sw = arm << quarter
    sw.ymin = +wwg/2+coupler_margin
    sw.xmax = 0

    se = arm << quarter
    se = se.mirror()
    se.connect("o2", sw.ports["o2"])

    arm.add_port("o1", port=sw.ports["o1"])
    arm.add_port("o2", port=se.ports["o1"])
    south = c << arm

    ubend1 = c << ubend
    ubend2 = c << ubend
    ubend1.connect("o1", south.ports["o1"])
    ubend2.connect("o2", south.ports["o2"])

    north = (c << arm)
    # north.mirror(p1=(0,0),p2=(1,0))
    north.connect("o2", ubend1.ports["o2"])

    xmin = c.xmin
    ymin = c.ymin
    ymax = c.ymax

    leg1 = c << gf.components.straight(length=l, width=wwg, layer=core_layer)
    # c.add_port(name="o1", center=(0, 0), width=wwg, orientation=0)
    leg1.ymax = ymin-coupler_margin
    leg1.xmin = xmin

    leg2 = c << gf.components.straight(length=l, width=wwg, layer=core_layer)
    leg2.ymin = ymax+coupler_margin
    leg2.xmin = xmin

    c.add_label("zigzag resonator", position=c.center)
    c.show()
    return c
