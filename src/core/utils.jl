function getdimsperm(frame)
    r = zeros(Int, 3)
    vs = vcat(collect(eachcol(Matrix(I, 3, 3))), collect(eachcol(-Matrix(I, 3, 3))))
    for i = 1:3
        j = findfirst(vs) do v
            all(isapprox.(frame[:, i], v; atol=0.01))
        end
        isnothing(j) && return nothing
        if j > 3
            j = 3 - j
        end
        r[i] = j
    end
    r
end


function get_polarization(u)
    # u=flatten(u)
    if haskey(u, :E)
        if length(u.E) == 2
            return :TE
        elseif length(u.H) == 2
            return :TM
        end
    else
        if haskey(u, :Hy)
            return :TE
        elseif haskey(u, :Ey)
            return :TM
        end
    end
    nothing
end

function add_current_keys!(d::AbstractDict)
    for k in keys(d)
        if startswith(String(k), "E")
            s = String(k)[end]
            # d[Symbol("J$s")] = d[Symbol("P$s")] = d[k]
            d[Symbol("J$s")] = d[k] |> deepcopy
            # d[Symbol("P$s")] = d[k] |> deepcopy
        end
    end
    d
end

function groupkeys(d)
    r = OrderedDict{Symbol,Any}()
    for k in keys(d)
        pre = String(k)[1] |> Symbol
        if !haskey(r, pre)
            r[pre] = OrderedDict{Symbol,Any}()
        end
        r[pre][k] = d[k]
    end
    OrderedDict{Any,Any}([k => NamedTuple(pairs(v)) for (k, v) = pairs(r)])
    # r
end

_make_field_deltas(d::Real, a...) = d
function _make_field_deltas(d, N, field_boundvals, sizes, i, isdiff=false)
    NamedTuple([k => begin
        if isdiff * isnothing(v[i, 1])
            d = vcat(d[1], (d[1:end-1] .+ d[2:end]) ./ 2)
        end
        sel = i .== 1:N
        reshape(d, Tuple(1 - sel + sizes[k][i] * sel))
    end for (k, v) = pairs(field_boundvals)])
end
function cloud(v)
    groups = []
    group = []
    for (i, x) = enumerate(v)
        if isempty(group) || x / v[group[1]] < 1.2
            push!(group, i)
        else
            push!(groups, group)
            group = [i]
        end
    end
    if !isempty(group)
        push!(groups, group)
    end
    groups
end
quasimax(x::Number, l) = abs(x)
function quasimax(a, l=0.8)
    a = abs.(a)
    a = filter(>(0), a .* (a .> mean(a)))
    isempty(a) && return 0
    quantile(a, l)
end

const SEP = "|"
_writejsonnp(_, x) = x
function _writejsonnp(pre, d::Union{AbstractDict,NamedTuple},)
    kvmap(d) do k, v
        k => _writejsonnp("$(pre)$SEP$k", v)
    end
end
function _writejsonnp(pre, a::AbstractArray)
    map(enumerate(a)) do (i, a)
        _writejsonnp("$(pre)$SEP$i", a)
    end
end
function _writejsonnp(pre, a::AbstractArray{<:Number},)
    name = "$(pre).npy"
    npzwrite(name, a)
    basename(name)
end
function writejsonnp(path, d::Union{AbstractDict,NamedTuple})
    DATA = randstring()
    mkpath(DATA)

    open(path, "w") do io
        write(io, JSON.json(kvmap(d) do k, v
                k => _writejsonnp(string(joinpath(DATA, string(k))), v)
            end, 4))
    end

    TAR = "$DATA.tar"
    Tar.create(DATA, TAR)

    _TAR = joinpath(dirname(path), ".$(basename(path)).tar")
    cp(TAR, _TAR, force=true)

    rm(TAR; force=true)
    rm(DATA; recursive=true)
end

_readjsonnp(pre, x::Union{AbstractDict,NamedTuple}) = OrderedDict([Symbol(k) => _readjsonnp("$pre$SEP$k", v) for (k, v) in pairs(x)])
_readjsonnp(pre, x::AbstractVector) = [_readjsonnp("$pre$SEP$i", x) for (i, x) in enumerate(x)]
function _readjsonnp(pre, x)
    endswith(string(x), ".npy") && return npzread("$pre.npy")
    x
end


function readjsonnp(path)
    _TAR = joinpath(dirname(path), ".$(basename(path)).tar")
    TAR = ".$(basename(path)).tar"
    cp(_TAR, TAR; force=true)

    DATA = randstring()
    mkpath(DATA)
    Tar.extract(TAR, DATA)

    r = OrderedDict([
        Symbol(k) => _readjsonnp(
            string(joinpath(DATA, k)), v
        )
        for (k, v) in pairs(JSON.parse(read(open(path), String); dicttype=OrderedDict))
    ])

    rm(DATA; recursive=true)
    rm(TAR; force=true)

    r
end