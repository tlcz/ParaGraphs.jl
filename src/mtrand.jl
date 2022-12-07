using  Random, Base.Threads
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

A multithreaded version of [`Random.randperm!`](@ref)

Construct in `A` a random permutation of length `length(A)`.
Arg `mask` determines number of parallel partitions to be used.
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
function mt_randperm!(r::TaskLocalRNG, v::Vector{T}, mask::Tu) where {T<:Integer, Tu<:Union{UInt8, UInt16}}
    nparts = mask + 1
    @assert ispow2(nparts) "invalid mask $(mask)"

    n  = length(v)
    s  = Random.SamplerType{Tu}()
    r0 = copy(r)

    counts = zeros(Int, nparts + 64, nthreads())

    # Random.seed!(r, seed)
    @threads :static for i in 1:n
        local tid, pid = threadid(), rand(r, s) & mask + 1
        @inbounds counts[pid, tid] += 1
    end

    prev = 0
    for p = 1:nparts, tid = 1:nthreads()
        @inbounds prev = counts[p, tid] += prev
    end

    copy!(r, r0)
    @threads :static for i in 1:n
        local tid, pid = threadid(), rand(r, s) & mask + 1
        @inbounds local ix = counts[pid, tid]
        @inbounds v[ix] = i
        @inbounds counts[pid, tid] -= 1
    end

    counts[nparts+1, 1] = n
    @threads :static for pid in 1:nparts
        @inbounds local chunk = view(v, counts[pid,1]+1:counts[pid+1,1])
        _shuffle!(r, chunk)
    end
    v
end

function mt_randperm!(r::TaskLocalRNG, v::Vector{T}) where {T<:Integer}
    nparts = (length(v) * sizeof(T)) >> 21
    nparts == 0 && return randperm!(v)
    nparts = nextpow(2, nparts)
    mask = nparts <= 256 ? UInt8(nparts-1) : UInt16(nparts-1)
    mt_randperm!(r, v, mask)
end

mt_randperm!(v::Vector{T}, mask::Tu) where {T<:Integer, Tu<:Union{UInt8, UInt16}} = mt_randperm!(default_rng(), v, mask)
mt_randperm!(v::Vector{T}) where {T<:Integer} = mt_randperm!(default_rng(), v)
