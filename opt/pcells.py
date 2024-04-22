import gdsfactory as gf


def phase_shifter(l, wwg, nbends, ubend):
    c = gf.Component()
    s = c << gf.path.extrude(gf.path.straight(length=l), width=wwg)
    for i in range(1, nbends):
        u_ = (c << ubend)
        p1, p2 = i % 2 == 1 ? "o1": "o2", i % 2 == 1 ? "o2": "o1"
        u_.connect(p1, s.ports["o2"])
        s_ = (c << gf.path.extrude(gf.path.straight(length=l), width=wwg))
        s_.connect("o1", u_.ports[p2])
        s = s_
    c.show()
    return c
