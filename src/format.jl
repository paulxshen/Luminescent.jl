
function shorten_key(k)
    s = split(k, ",")
    p1, p2 = port_number.(s)
    mn1, mn2 = mode_number.(s)
    if mn1 == mn2 == 0
        return "$p1,$p2"
    end
    return nothing
end
function simplify_sparams(s)
    k = sort(keys(s))
    k = k[round((length(k) + 1) / 2)]
    dict([shorten_key(k) => v for (k, v) = pairs(s[k]) if !isnothing(shorten_key(k))])
end
function sparam_family(sparams)
    # sparams = simplify_sparams(detailed_sparams)
    tparams = Porcupine.apply(abs2, sparams)
    # detailed_tparams = Porcupine.apply(abs2, detailed_sparams)
    # sparam_phasors = Porcupine.apply(sparams) do z
    #     (mag=abs(z), phase=angle(z))
    # end
    # println(tparams)
    # sparams = apply(reim, sparams)
    sol = (; sparams, tparams,)# detailed_sparams, detailed_tparams)#sparam_phasors)
end
function port_number(port)
    s = split(string(port), "@")[1]
    if s[1] == 'o'
        s = s[2:end]
    end
    return parse(Int, s)
end

function mode_number(port)
    l = split(string(port), "@")
    return if length(l) == 1
        0
    else
        parse(Int, l[2])
    end
end


# function long_sparam_key(s):
#     o, i = s.split(",")
#     po, pi = port_number(o), port_number(i)
#     mo, mi = mode_number(o), mode_number(i)
#     return f"o{po}@{mo},o{pi}@{mi}"
