function load(s, _dx)
    @load "$s.bson" base design_start design_sz dx
    L = size(base) .* dx
    ratio = dx / _dx
    l = (58.0f0 / 108) * L[1]
    design_sz = round.(Int, design_sz .* ratio)
    design_start = reindex(design_start, ratio)
    base = imresize(base; ratio)
    ports = ([0, l], [l, 0])
    (; base, design_start, design_sz, ports)
end

index(v, dx) = round.(Int, v ./ dx .+ 1)
reindex(i, ratio) = round.(Int, (i .- 1) .* ratio .+ 1)
Base.size(x::NamedTuple) = (length(x),)