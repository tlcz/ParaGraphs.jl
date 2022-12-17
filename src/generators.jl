using SparseArrays
import SparseArrays: AbstractSparseMatrixCSC, getcolptr

function config_model(degs::Vector{T}) where T<:Integer
    n, s = length(degs), sum(degs)
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
    @time perm = mt_randperm!(similar(rows))

    @info "generate matching"
    @time matching(rows, perm)
    @info "creating graph"
    @time m = SparseMatrixCSC(n, n, colptr, rows, ones(Int8, s));
    @info "sorting columns"
    @time sortCols3!(m)
    @info "trimming zeros"
    @time dropzeros!(m)
end

function matching(rows::Vector{T}, perm::Vector{T}) where T<:Integer
    @inbounds @threads for i in 1:2:length(perm)
        # p = perm[i] + perm[i+1]
        @inbounds local p, q = perm[i], perm[i+1]
        @inbounds rows[p], rows[q] = rows[q], rows[p]
    end
end

function _mt_populate!(dest::AbstractArray{T}, colptr::AbstractArray{T}) where {T<:Integer}
    @threads :static for i in 1:length(colptr)-1
        @inbounds local s,e = colptr[i],colptr[i+1]-1
        @inbounds dest[s:e] .= i
    end
end

function sortCols3!(A::AbstractSparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
    m, n = size(A)
    # colptr = getcolptr(A);
    rowval = rowvals(A);
    nzval = nonzeros(A)

    @inbounds @threads :static for i = 1:n
        local nzr = nzrange(A, i)
        local numrows = length(nzr)

        numrows <= 1 && continue

        alg = numrows <= 16 ? Base.Sort.InsertionSort : Base.Sort.QuickSort
        @inbounds sort!(view(rowval, nzr), 1, length(nzr), alg, Base.Order.Forward)

        # coalesce duplicates
        local rprev,ix = 0,0
        @inbounds for j = nzr
            (rowval[j]>rprev) && (ix=j)
            rowval[j]==rprev && (nzval[ix]+=nzval[j]; nzval[j]-=nzval[j])
            rprev = rowval[j]
        end
    end

    return A
end
