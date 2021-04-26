module ROIViews

export ROIView, expand_size

 # T refers to the result type. N to the dimensions of the final array, and M to the dimensions of the raw array
struct ROIView{T, N, M, Z, AA<:AbstractArray} <: AbstractArray{T, N}  # this type has one dimension more than input
    # stores the data. 
    parent::AA  # the array from which the ROIs are referenced
    # output size of the array 
    views::NTuple{Z, SubArray{T, M, AA, NTuple{M, UnitRange{Int}}, false}}  # a tuple of views indexed by the last index
    num_ROIs::Int
    pad_val::T
    mysize::NTuple{N,Int}

    # Constructor function
    function ROIView{T, N, M, Z}(data::AA, views::NTuple{Z,SubArray}, pad_val::T) where {T,N,M,Z,AA}
        sz = (expand_size(size(views[1]), size(data)) ..., Z)
        return new{T, N, M, Z, AA}(data, views, Z, pad_val, sz) 
    end
end


"""
    ROIView([T], data::F, tile_size::NTuple{N,Int}, tile_overlap::NTuple{N,Int}) where {N,F}

Creates an M+1 dimensional view of an array by stacking regions of interest (ROIs) along the first unused dimension. 

`data`. the inputdata hosting the ROIs.

`ROI_size`. The common size of all ROIs.

# Examples
```jldoctest
```
"""
function ROIView(data::AA, center_pos, ROI_size::NTuple; pad_val=0) where {AA}
    of_starts = Tuple(Tuple(pos .- (ROI_size .÷ 2)) for pos in center_pos)
    DataDims = ndims(data)
    of_starts = Tuple(expand_zeros(s, Val(DataDims)) for s in of_starts)
    ROI_size = Tuple(expand_size(ROI_size, size(data)))
    views = Tuple(view(data,(ofs[d]+1:ofs[d]+ROI_size[d] for d=1:DataDims)...) for ofs in of_starts)
    return ROIView{eltype(data), DataDims+1, DataDims, length(of_starts)}(data, views, convert(eltype(data),pad_val)) 
end

function ROIView(data::AA, center_pos::Matrix, ROI_size::NTuple; pad_val=0) where {AA}
    return ROIView(data, eachslice(center_pos, dims=2), ROI_size, pad_val=pad_val) 
end

function expand_size(sz,sz2)
    dims1 = length(sz)
    dims2 = length(sz2)
    ((d<=dims1) ? sz[d] : sz2[d] for d in 1:dims2)
end

function expand_zeros(t, ::Val{N}) where N  # adds t1 to t2 as a tuple and returns t2[n] for n > length(t1)
    ntuple(x-> x ≤ length(t) ? t[x] : 0,Val(N))
end

# define AbstractArray function to allow to treat the generator as an array
# See https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array
function Base.size(A::ROIView)
    return A.mysize  # appends num_ROIs to the ROI_size
end

Base.similar(A::ROIView, ::Type{T}, size::Dims) where {T} = ROIView{T, ndims(A), A.num_ROIs}(A.parent, A.ranges, A.ROIsize)

# %24 = Base.getproperty(A, :parent)::AbstractMatrix{Float64}

# calculate the entry according to the index
# Base.getindex(A::IndexFunArray{T,N}, I::Vararg{B, N}) where {T,N, B} = return A.generator(I)

function expand_add(t1,t2)  # adds t1 to t2 as a tuple and returns t2[n] for n > length(t1)
    ((t+w for (t,w) in zip(t1,t2))..., (w for w in t2[length(t1)+1:end])...)
    # ((t+w for (t,w) in zip(t1,t2))...,t2[length(t1)+1:end]...)
end

# calculate the entry according to the index
function Base.getindex(A::ROIView{T,N}, I::Vararg{Int, N}) where {T,N}
    @boundscheck checkbounds(A, I...) # check acces for this view
    myview = A.views[last(I)]
    pos = Base.front(I)
    if Base.checkbounds(Bool, myview, pos...) # check acces in original array
        return Base.getindex(myview, pos... ) # pos
    else
        return A.pad_val
    end
end

# not supported
Base.setindex!(A::ROIView{T,N}, v, I::Vararg{Int,N}) where {T,N} = begin 
    @boundscheck checkbounds(A, I...)
    myview = A.views[last(I)]
    pos = Base.front(I)
    if Base.checkbounds(Bool, myview, pos...) # check acces in original array
        return Base.setindex!(myview, v, Base.front(I)... ) # pos
    else
        return A.pad_val
    end
end


end # module
