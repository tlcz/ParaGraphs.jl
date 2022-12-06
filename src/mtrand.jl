using  Random, Base.Threads
import Random: Sampler, default_rng, seed!

"""
    mt_randperm!([rng=GLOBAL_RNG,] A::Array{<:Integer}, mask::Union{UInt8, UInt16})

multithreaded version of [`randperm!`](@ref)

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
        shuffle!(r, chunk)
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
