struct Cat{T,N}<:AbstractArray{T,N}
    v
    dims
ends
sz
len
function Cat(a...;dims)
ends=cumsum(size.(a,dims))
sz=size(a[1])
sz=(sz[1:dims-1],ends[end],sz[dims+1:end])
len=sum(length.(a))
    new{eltype(a[1]),ndims(a[1])}(a,dims,ends,sz,len)
end
end
lcat(a...;dims)=Cat(a...;dims)
lvcat(a...)=lcat(a...;dims=1)
lhcat(a...)=lcat(a...;dims=2)
Base.size(m::Cat)=m.sz
Base.length(m::Cat)=m.len
Base.collect(m::Cat)=cat(m.v...;dims=m.dims)
Base.vec(m::Cat)=vcat(vec.(m.v)...;)
function overlap(ends,i::Int)
    vi=ss(ends,i)
    if vi>1
  return      [vi],[i-ends[vi-1]]
    end
     [vi],[i]
        end
function overlap(ends,r::)
    (vi1,),(i1,)=overlap(ends,r[1])
    (vi2,),(i2,)=overlap(ends,r[end])
    vi1:vi2,[if vi==vi1
    i1:ends[vi1]
elseif vi==vi2
    1:i2
else
    (:)
end for vi=vi1:vi2]
        end
function Base.getindex(m::Cat,I...)
@unpack v,dims,ends=m
vi,vr=overlap(ends,I[dims])
I1=I[1:dims-1]
I2=I[dims+1:end]
cat([v[i][I1...,r,I2...] for (i,r)=zip(vi,vr)]...;dims)
end