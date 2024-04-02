using VectorModesolver
using BSON: @save, @load

# function (ε::εtype)(x, y, w, h, wm, hm)
function main(dx, λ, w, h, wm, hm, ϵ1, ϵ2, neigs=1)
    tol = 1.0e-4
    # ex = Meta.parse("""
    #       function (ε::εtype)(x, y,)
    #           # lm=()
    #           if ($wm-$tol <= x < $w - $wm-$tol) && ($hm-$tol <= y < $h - $hm-$tol)
    #               return ($ϵ2, 0.0, 0.0, $ϵ2, $ϵ2)
    #           end

    #           return ($ϵ1, 0.0, 0.0, $ϵ1, $ϵ1)
    #       end
    #       """)
    # eval(ex)
    # sleep(2)
    ε = εtype()
    x = [i for i in 0:dx:w-tol]
    y = [i for i in 0:dx:h-tol]
    # tol = 1e-8
    boundary = (0, 0, 0, 0)
    # eval()
    solver = VectorialModesolver(λ, x, y, boundary, ε)
    modes = solve(solver, neigs, tol)

end
function rectmodes(dx, λ, wwg, hwg, ϵ1, ϵ2, wm=0.25, hm=0.25,)
    # wm = hsub = hclad = 0.255
    w = wwg + 2wm
    h = hwg + 2hm
    lb = [-w / 2, -hm]
    ub = [w / 2, hwg + hm]

    modes_ = main(dx, λ, w, h, wm, hm, ϵ1, ϵ2)
    modes = [NamedTuple([k => transpose(getfield(m, k),) for k = [:neff, :Ex, :Ey, :Ez, :Hx, :Hy, :Hz]]) for m = modes_]
    # @save fn modes λ dx ub lb hsub hwg wwg hclad ϵsub ϵclad ϵcore h w
    display(plot_mode_fields(modes_[1])
    )
    (; modes, λ, dx, ub, lb, hsub, hwg, wwg, hclad, ϵsub, ϵclad, ϵcore, h, w)
end
function linemodes(dx, λ, wwg, ϵ1, ϵ2, wm=0.25, hm=wm)
    # wm = hsub = hclad = 0.255
    w = wwg + 2wm
    lb = [-w / 2,]
    ub = [w / 2,]
    h = 4wwg + 2hm
    modes_ = main(dx, λ, w, h, wm, hm, ϵ1, ϵ2, 3)
    modes = [NamedTuple([k => transpose(getfield(m, k),) for k = [:neff, :Ex, :Ey, :Ez, :Hx, :Hy, :Hz]]) for m = modes_]
    modes = [
        begin
            @unpack Ex, Ey, Ez, Hx, Hy, Hz = m
            Ey, Ex, Hz = map([Ex, Ez, Hy]) do a
                transpose(sum(a, dims=2))
            end
            (; Ex, Ey, Hz)
        end for m = modes
    ]
    # @save fn modes λ dx ub lb hsub hwg wwg hclad ϵsub ϵclad ϵcore h w
    display(plot_mode_fields(modes_[2]))
    display(plot_mode_fields(modes_[3]))
    # display(plot_mode_fields(modes_[4]))
    # display(plot_mode_fields(modes_[5]))
    (; modes=[modes[3]], λ, dx, ub, lb, wwg, ϵclad, ϵcore, w)
end