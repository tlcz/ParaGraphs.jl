using SparseArrays
import Random: mt_randperm
# import SparseArrays: AbstractSparseMatrixCSC, getcolptr

function config_model(degs::Vector{T}) where T<:Integer
    n = length(degs)
    @assert n <= typemax(T)
    s = sum(degs)
    @assert iseven(s) "sum of degrees is odd"
    mn, mx = extrema(degs)
    @assert mn > 0 && mx <= n

    @info "preparing half-edges"
    # generate colptr
    colptr = similar(degs, n + 1)
    colptr[1] = 1; copyto!(colptr, 2, degs)
    cumsum!(colptr, colptr)
    @assert @inbounds colptr[end] == s + 1
    rows = similar(degs, s)
    @time _mt_populate!(rows, colptr)

    @info "generate permutation"
    @time perm = mt_randperm(s%T)
    @info "generate matching"
    @time matching(rows, perm) #, 100000000%T)

    @info "creating graph"
    @time m = SparseMatrixCSC(n, n, colptr, rows, ones(Int8, s));
    @info "sorting columns"
    @time sortCols3!(m)
    @info "trimming zeros"
    @time dropzeros!(m)
end

function halfperm(N::T) where T<:Integer
    @time b = bitrand(N)
    @time p = findall(b)
    @time l = length(p)
    @time resize!(p, N)
    @time map!(~, b, b)
    @time mt_shuffle!(view(p,l+1:N), findall(b))
    p
end

function matching(rows::Vector{T}, perm::Vector{T}) where T<:Integer
    perm = reshape(perm, :, 2)
    @inbounds @threads for i in 1:size(perm,1)
        # p = perm[i] + perm[i+1]
        @inbounds local p, q = perm[i,1], perm[i,2]
        @inbounds rows[p], rows[q] = rows[q], rows[p]
    end
end

function matching(rows::Vector{T}, perm::Vector{T}, k::T) where T<:Integer
    perm = reshape(perm, :, 2)
    r = [x:x+k-1 for x in 1:k:length(rows)]

    @threads :static for pr in r for qr in r
    for i in 1:size(perm,1)
        @inbounds local p, q = perm[i,1], perm[i,2]
        p ∈ pr && q ∈ qr && @inbounds rows[p], rows[q] = rows[q], rows[p]
    end end end
end

function matching(rows::Vector{T}, colptr::Vector{T}, perm::Vector{T}) where T<:Integer
    perm = reshape(perm, :, 2)
    @inbounds @threads for i in 1:size(perm,1)
        # p = perm[i] + perm[i+1]
        @inbounds local p, q = perm[i,1], perm[i,2]
        @inbounds rows[p], rows[q] = searchsortedlast(colptr, q), searchsortedlast(colptr, p)
    end
end

function _mt_populate!(dest::AbstractArray{T}, colptr::AbstractArray{T}) where {T<:Integer}
    @threads :static for i in 1:length(colptr)-1
        @inbounds local s,e = colptr[i],colptr[i+1]-1
        @inbounds dest[s:e] .= i
    end
end

function sortCols3!(A::AbstractSparseMatrix{Tv,Ti}) where {Tv,Ti}
    m, n = size(A)
    # colptr = getcolptr(A);
    rowval = rowvals(A);
    nzval = nonzeros(A)

    @inbounds @threads :static for i = 1:n
        local nzr = nzrange(A, i)
        local nrows = length(nzr)

        nrows <= 1 && continue

        alg = nrows <= 16 ? Base.Sort.InsertionSort : Base.Sort.QuickSort
        @inbounds sort!(view(rowval, nzr), 1, length(nzr), alg, Base.Order.Forward)

        # coalesce duplicates
        local rprev, ix = 0%Tv, 0
        @inbounds for j = nzr
            @inbounds (rowval[j]>rprev) && (ix=j)
            @inbounds rowval[j]==rprev && (nzval[ix]+=nzval[j]; nzval[j]-=nzval[j])
            @inbounds rprev = rowval[j]
        end
    end

    return A
end

function dupa(v)
   n = length(v)%eltype(v)
   b = Base.OneTo(n)
   for i in b
       ix = rand(b)
       x = v[i]
       v[ix] = i
   end
end
