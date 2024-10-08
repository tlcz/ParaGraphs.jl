using  Random, Base.Threads
import Base: BitInteger
import Random: Sampler, default_rng, seed!
import Random: require_one_based_indexing, ltm52, shuffle!

"""
    _shuffle(r::AbstractRNG, a::AbstractArray)

A sligthly faster version of [`Random.shuffle!`](@ref)
"""
function _shuffle!(r::AbstractRNG, a::AbstractArray)
    require_one_based_indexing(a)
    n = length(a)
    n <= 1 && return a # nextpow below won't work with n == 0
    @assert n <= Int64(2)^52
    mask = nextpow(2, n) - 1
    for i = n:-1:2
        (mask >> 1) == i && (mask >>= 1)
        j = 1 + rand(r, ltm52(i, mask))
        @inbounds a[i], a[j] = a[j], a[i]
    end
    return a
end

"""
    mt_randperm!([rng=GLOBAL_RNG,] A::Array{<:Integer}, mask<:Union{UInt8, UInt16})

A multithreaded version of [`Random.randperm!`](@ref).
It works by first partitioning input uniformly at random
and then shuffling resulting partitions in parallel.

Construct in `A` a random permutation of length `length(A)`.
Arg `mask` determines number of parallel partitions (mask + 1) to be used.
Optional arg `rng` specifies a random number generator (see [Random Numbers](@ref)).

# Examples
```jldoctest
julia> r = TaskLocalRNG();

julia> Random.seed!(r, 1234);

julia> mt_randperm!(r, Vector{Int}(undef, 8), 0x3)
8-element Vector{Int64}:
 2
 1
 7
 4
 3
 5
 6
 8
```
"""
function mt_randperm!(r::TaskLocalRNG, A::Array{T}, mask::Tu) where {T<:BitInteger, Tu<:Union{UInt8, UInt16}}
    # determine number of partitions
    nparts = mask + 1
    @assert ispow2(nparts) "invalid mask $(mask)"

    n  = length(A)
    s  = Random.SamplerType{Tu}()

    # array storing partition boundaries (partition by thread)
    # partitions map to rows (pid)
    # threads map to cols (tid)
    # we add some extra space for each thread
    # to avoid inter-core synchronization
    # of cached adjacent mem addresses
    counts = zeros(Int, nparts + 64, nthreads())

    # save initial random state
    r0 = copy(r)
    # 1st pass
    # assign input to partitions uniformly at random
    # calculate space claimed by each thread in every partition
    @threads :static for i in 1:n
        local tid, pid = threadid(), rand(r, s) & mask + 1
        @inbounds counts[pid, tid] += 1
    end

    # cumsum partition sizes
    # to mark boundaries of space claimed by each thread in every partition
    # note that the 1st column will contain boundaries of entire partitions
    prev = 0
    for pid = 1:nparts, tid = 1:nthreads()
        @inbounds prev = counts[pid, tid] += prev
    end
    # mark the end of the last partition
    counts[nparts + 1, 1] = n

    # recover random state
    copy!(r, r0)
    # 2nd pass
    # scatter input accross partitions uniformly at random
    # note that input distribution is identical as in the 1st pass
    # since we recovered the initial random state
    @threads :static for i in 1:n
        local tid, pid = threadid(), rand(r, s) & mask + 1
        @inbounds A[counts[pid, tid]] = i
        @inbounds counts[pid, tid] -= 1
    end

    # input is partitioned
    # shuffle partitions in parallel
    @threads :static for pid in 1:nparts
        @inbounds local chunk = view(A, counts[pid, 1] + 1:counts[pid + 1, 1])
        _shuffle!(r, chunk)
    end
    A
end

function mt_randperm!(r::TaskLocalRNG, A::Array{T}) where {T<:Integer}
    nthreads() == 1 && return randperm!(A)
    nparts = (length(A) * sizeof(T)) >> 21
    nparts == 0 && return randperm!(A)
    nparts = nextpow(2, nparts)
    mask = nparts <= typemax(UInt8) + 1 ? UInt8(nparts - 1) : UInt16(nparts - 1)
    mt_randperm!(r, A, mask)
end

mt_randperm!(A::Array{T}, mask::Tu) where {T<:Integer, Tu<:Union{UInt8, UInt16}} = mt_randperm!(default_rng(), A, mask)
mt_randperm!(A::Array{T}) where {T<:Integer} = mt_randperm!(default_rng(), A)
