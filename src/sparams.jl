function vec2smatrix(v)
    n = length(v)
    map(CartesianIndices((n, n))) do I
        I = Tuple(I)
        i, j = I
        if j == 1
            v[i]
        elseif i == 1
            v[j]
        else
            1
        end
    end
end

function sparam_family(S)
    # S = simplify_sparams(detailed_sparams)
    T = fmap(abs2, S)
    dB = fmap(T) do x
        10log10(x)
    end
    # detailed_tparams = Porcupine.apply(abs2, detailed_sparams)
    phasors = fmap(S) do z
        (mag=abs(z), phase=rad2deg(angle(z)))
    end
    sol = (; S, T, phasors, dB) |> pairs |> OrderedDict
end

