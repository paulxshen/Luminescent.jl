
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
    k = k[round(Int, (length(k) + 1) / 2)]
    dict([shorten_key(k) => v for (k, v) = pairs(s[k]) if !isnothing(shorten_key(k))])
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
