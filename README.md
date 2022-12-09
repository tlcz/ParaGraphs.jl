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

Generate config model graph - a toy example
```julia
julia> using Random, LinearAlgebra, ParaGraphs

julia> Random.seed!(123);

julia> d = Int32.([8,7,6,5,4,3,2,1]);

julia> G = config_model(d);
[ Info: generate half-edges
  0.000000 seconds
[ Info: generate permutation
  0.000001 seconds (1 allocation: 208 bytes)
[ Info: generate matching
  0.000107 seconds (50 allocations: 4.391 KiB)
[ Info: creating graph
  0.000001 seconds (1 allocation: 96 bytes)
[ Info: sorting
  0.000065 seconds (50 allocations: 4.641 KiB)
[ Info: cleansing
  0.000001 seconds
8×8 SparseMatrixCSC{Int8, Int32} with 24 stored entries:
 ⋅  3  ⋅  2  2  1  ⋅  ⋅
 3  ⋅  2  1  1  ⋅  ⋅  ⋅
 ⋅  2  ⋅  1  ⋅  2  ⋅  1
 2  1  1  ⋅  ⋅  ⋅  1  ⋅
 2  1  ⋅  ⋅  ⋅  ⋅  1  ⋅
 1  ⋅  2  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  1  1  ⋅  ⋅  ⋅
 ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅

julia> issymmetric(G)
true

julia> sum(G, dims=1)[:] == d
true

julia> G[:,2]
8-element SparseVector{Int8, Int32} with 4 stored entries:
  [1]  =  3
  [3]  =  2
  [4]  =  1
  [5]  =  1

julia> findnz(G)
(Int32[1, 2, 3, 4, 5, 6, 8, 1, 3, 4  …  1, 2, 3, 4, 1, 3, 7, 4, 6, 1], Int32[1, 1, 1, 1, 1, 1, 1, 2, 2, 2  …  5, 5, 5, 5, 6, 6, 6, 7, 7, 8], Int8[2, 1, 1, 1, 1, 1, 1, 1, 3, 2  …  1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
```
extract multi-edges
```julia
julia> G .> 1
8×8 SparseMatrixCSC{Bool, Int32} with 10 stored entries:
 ⋅  1  ⋅  1  1  ⋅  ⋅  ⋅
 1  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  1  ⋅  ⋅  ⋅  1  ⋅  ⋅
 1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
```

a serious example
```julia
julia> n = Int(1e7);

julia> Random.seed!(321);

julia> d = rand(Int32(1):Int32(400), n);

julia> @time G = config_model(d);
[ Info: generate half-edges
  2.805701 seconds
[ Info: generate permutation
 11.066151 seconds (226 allocations: 7.470 GiB, 0.18% gc time)
[ Info: generate matching
 17.654372 seconds (48 allocations: 4.328 KiB)
[ Info: creating graph
  0.809052 seconds (2 allocations: 1.867 GiB, 8.01% gc time)
[ Info: sorting
  9.066615 seconds (49 allocations: 4.609 KiB)
[ Info: cleansing
  1.599250 seconds
 43.037295 seconds (856 allocations: 16.845 GiB, 0.21% gc time)

julia> sum(G, dims=1)[:] == d
true

julia> varinfo()
 name                    size summary                                                                                   
 –––––––––––––––– ––––––––––– ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
 d                 38.147 MiB 10000000-element Vector{Int32}                                                            
 G                  9.375 GiB 10000000×10000000 SparseMatrixCSC{Int8, Int32} with 2005176310 stored entries
