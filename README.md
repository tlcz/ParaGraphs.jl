# ParaGraphs.jl
Parallel random graph generators for Julia

Generate config model graph - a toy example (8 threads)
```julia
julia> using Random, LinearAlgebra, ParaGraphs

julia> Random.seed!(123);

julia> d = Int32.([8,7,6,5,4,3,2,1]);

julia> @time G = config_model(d)
[ Info: preparing half-edges
  0.000088 seconds (50 allocations: 4.266 KiB)
[ Info: generate permutation
  0.000001 seconds (1 allocation: 208 bytes)
[ Info: generate matching
  0.000105 seconds (49 allocations: 4.359 KiB)
[ Info: creating graph
  0.000001 seconds (1 allocation: 96 bytes)
[ Info: sorting columns
  0.000116 seconds (49 allocations: 4.609 KiB)
[ Info: trimming zeros
  0.000001 seconds
  0.000672 seconds (687 allocations: 41.336 KiB)
8×8 SparseArrays.SparseMatrixCSC{Int8, Int32} with 26 stored entries:
 2  2  2  1  ⋅  1  ⋅  ⋅
 2  ⋅  2  ⋅  1  1  1  ⋅
 2  2  ⋅  2  ⋅  ⋅  ⋅  ⋅
 1  ⋅  2  2  ⋅  ⋅  ⋅  ⋅
 ⋅  1  ⋅  ⋅  ⋅  1  1  1
 1  1  ⋅  ⋅  1  ⋅  ⋅  ⋅
 ⋅  1  ⋅  ⋅  1  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅

julia> issymmetric(G)
true

julia> sum(G, dims=1)[:] == d
true

julia> G[:,2]
8-element SparseArrays.SparseVector{Int8, Int32} with 5 stored entries:
  [1]  =  2
  [3]  =  2
  [5]  =  1
  [6]  =  1
  [7]  =  1

julia> using SparseArrays
julia> findnz(G)
(Int32[1, 2, 3, 4, 6, 1, 3, 5, 6, 7  …  2, 6, 7, 8, 1, 2, 5, 2, 5, 5], Int32[1, 1, 1, 1, 1, 2, 2, 2, 2, 2  …  5, 5, 5, 5, 6, 6, 6, 7, 7, 8], Int8[2, 2, 2, 1, 1, 2, 2, 1, 1, 1  …  1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
```
extract multi-edges
```julia
julia> G .> 1
8×8 SparseMatrixCSC{Bool, Int32} with 10 stored entries:
 1  1  1  ⋅  ⋅  ⋅  ⋅  ⋅
 1  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅
 1  1  ⋅  1  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  1  1  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
```

a serious example
```julia
julia> n = Int(1e7);

julia> Random.seed!(321);

julia> d = rand(Int32(1):Int32(400), n);

julia> @time G = config_model(d);
[ Info: preparing half-edges
  0.597114 seconds (51 allocations: 4.297 KiB)
[ Info: generate permutation
 11.506196 seconds (149 allocations: 7.470 GiB, 0.02% gc time)
[ Info: generate matching
 17.822963 seconds (51 allocations: 4.422 KiB)
[ Info: creating graph
  0.765683 seconds (2 allocations: 1.867 GiB, 1.24% gc time)
[ Info: sorting columns
  8.967472 seconds (53 allocations: 4.734 KiB)
[ Info: trimming zeros
  1.629879 seconds
 41.321492 seconds (846 allocations: 16.845 GiB, 0.03% gc time)

julia> sum(G, dims=1)[:] == d
true

julia> varinfo()
 name                    size summary                                                                                   
 –––––––––––––––– ––––––––––– ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
 d                 38.147 MiB 10000000-element Vector{Int32}                                                            
 G                  9.375 GiB 10000000×10000000 SparseMatrixCSC{Int8, Int32} with 2005176310 stored entries
```

a faster multi-threaded randperm! on board
```julia
julia> n = 1000^3;
julia> v = Vector{Int32}(undef, n);

julia> using ParaGraphs
julia> @time mt_randperm!(v);
  3.645866 seconds (77 allocations: 72.531 KiB)

julia> isperm(v)
true

julia> using Random
julia> @time randperm!(v);
 31.672325 seconds (14 allocations: 784 bytes, 0.01% compilation time)
```
