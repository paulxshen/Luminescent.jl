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