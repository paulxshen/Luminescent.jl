using VectorModesolver
using BSON: @save, @load

# function (ε::εtype)(x, y, w, h, wm, hm)
function main(dx, w, h, wm, hm)
    ex = Meta.parse("""
          function (ε::εtype)(x, y,)
              # lm=()
              if ($wm <= x <= $w - $wm) && ($hm <= y <= $h - $hm)
                  return ($ϵ2, 0.0, 0.0, $ϵ2, $ϵ2)
              end
              
              return ($ϵ1, 0.0, 0.0, $ϵ1, $ϵ1)
          end
          """)
    eval(ex)
    ε = εtype()
    x = [i for i in 0:dx:w]
    y = [i for i in 0:dx:h]
    neigs = 1
    tol = 1e-8
    boundary = (0, 0, 0, 0)
    # eval()
    solver = VectorialModesolver(λ, x, y, boundary, ε)
    modes = solve(solver, neigs, tol)

end
function savemodes(fn, dx, λ, wwg, hwg, wm, hm, ϵ1, ϵ2,)
    # wm = hsub = hclad = 0.25
    w = wwg + 2wm
    h = hwg + 2hm
    lb = [-w / 2, -hm]
    ub = [w / 2, hwg + hm]

    modes_ = main(dx, w, h, wm, hm)
    modes = [NamedTuple([k => transpose(getfield(m, k),) for k = [:neff, :Ex, :Ey, :Ez, :Hx, :Hy, :Hz]]) for m = modes_]
    # @save fn modes λ dx ub lb hsub hwg wwg hclad ϵsub ϵclad ϵcore h w
    display(plot_mode_fields(modes_[1])
    )
    (; modes, λ, dx, ub, lb, hsub, hwg, wwg, hclad, ϵsub, ϵclad, ϵcore, h, w)
end