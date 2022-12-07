# ParaGraphs.jl
Parallel graph algorithms for Julia
```julia
julia> using ParaGraphs

julia> n = 1000^3;

julia> v = Vector{Int32}(undef, n);

julia> @time mt_randperm!(v);
  3.645866 seconds (77 allocations: 72.531 KiB)

julia> isperm(v)
true
```
