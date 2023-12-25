using DelimitedFiles

function load(dir, s, _dx;)
    @load "$dir/$s" base design_start design_sz dx
    L = size(base) .* dx
    lp = (58.0f0 / 108) * L[1]
    ld = design_sz[1] * dx
    ratio = dx / _dx
    design_sz = round.(Int, design_sz .* ratio)
    design_start = reindex(design_start, ratio)
    base = imresize(base; ratio)
    ports = Any[[0, lp], [lp, 0]]
    wwg = 0.2
    lwg = 1.4
    TE0 = eachcol(readdlm("$dir/TE0.txt", ' ', Float32, '\n'))
    (; base, design_start, design_sz, ports, wwg, lwg, ld, lp, TE0)
end

index(v, dx) = round.(Int, v ./ dx .+ 1)
reindex(i, ratio) = round.(Int, (i .- 1) .* ratio .+ 1)
Base.size(x::NamedTuple) = (length(x),)