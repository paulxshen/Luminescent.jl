struct SubpixelAveraging
    v
end

function _apply_subpixel_averaging(s::SubpixelAveraging, a::AbstractArray)
    d = ndims(a)
    # for (i, (do_smooth, pad_left, pad_right)) in enumerate(s.v)
    for (i, (l, r)) = enumerate(eachcol(s.v))
        select = i .== (1:d) |> Tuple
        if l == -1
            a = pad(a, :replicate, select, 0)
        end
        if r == 1
            a = pad(a, :replicate, 0, select)
        end
        if l == -1 == r
            a = (2selectdim(a, i, 2:(size(a, i))) - diff(a, dims=i)) / 2
        elseif l == 1 == r
            a = (2selectdim(a, i, 1:(size(a, i)-1)) + diff(a, dims=i)) / 2
        end
    end
    a
end


function apply_subpixel_averaging(sas, gs)
    sas = ignore_derivatives() do
        sas
    end
    namedtuple([k => _apply_subpixel_averaging(sas[k], values(gs)[findfirst(keys(gs)) do gk
        startswith(string(k), string(gk))
    end]) for k = keys(sas)])
end


function _apply_geometry_padding(p::AbstractVector{<:OutPad}, a)
    for p = p
        @unpack l, r, b = p
        a = pad(a, b, l, r)
    end
    a
end

function apply_geometry_padding(gps, gs)
    namedtuple([k => _apply_geometry_padding(gps[k], gs[k]) for k = keys(gs)])
end

function _apply_field_padding(p::AbstractVector{<:InPad}, a::AbstractArray; nonzero_only=false)
    if nonzero_only
        p = filter(p) do p
            p.b != 0
        end
    end
    isempty(p) && return a

    if length(p) == 1 && !isnothing(p[1].m)
        return a .* p[1].m
    end

    a_ = Buffer(a)
    # a_ .= a
    a_[axes(a)...] = a
    for p = p
        @unpack l, r, b, m = p
        if isnothing(m)
            a_ = pad!(a_, b, l, r;)
        else
            a_[axes(a)...] = a .* m
        end
    end

    copy(a_)
end

function apply_field_padding(fps, fs; kw...)
    # namedtuple([k => apply(fps[k], fs[k]; kw...) for k = keys(fs)])
    dict([k => _apply_field_padding(fps[k], fs[k]; kw...) for k = keys(fs)])
end

function _mark(v::AbstractVector, a::AbstractArray)
    l = sum(v) do p
        p.l
    end
    r = sum(v) do p
        p.r
    end
    PaddedArray(a, l, r)
end

function mark(p, kw)
    namedtuple([k => _mark(p[k], kw[k]) for k = keys(kw)])
end

function unmark(kw)
    namedtuple([k => array(kw[k]) for k = keys(kw)])
end
