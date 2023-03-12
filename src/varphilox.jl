@inline PHILOX_M2x_0(::Type{UInt64}) = 0xD2B74407B1CE6E93
@inline PHILOX_M2x_0(::Type{UInt32}) = 0xd256d193

for (w, T, Td) in ((32, UInt32, UInt64), (64, UInt64, UInt128))
    @eval @inline function philox_mulhilo(a::$Td, b::$T)
        product = (a % $Td) * (b % $Td)
        hi = (product >> $w) % $T
        lo = product % $T
        (hi, lo)
    end
end

function varphilox(val::UInt64, k::Vector{UInt32}, lbits=8, rbits=8)
    lmask = ( one(UInt32) << lbits ) - one(UInt32)
    rmask = ( one(UInt32) << rbits ) - one(UInt32)
    state0, state1 = UInt32(val >> rbits), UInt32(val & rmask)
    for key in k
        hi, lo = philox_mulhilo( PHILOX_M2x_0(UInt64), state0) # TODO: UInt64
        lo = ( lo << ( rbits - lbits ) ) | state1 >> lbits
        state0 = ( ( hi ⊻ key ) ⊻ state1 ) & lmask
        state1 = lo & rmask
    end
    # Combine the left and right sides together to get result
    return (UInt64(state0) << rbits | UInt64(state1))
end

function dupa(N::Integer, R::Integer, bits::Integer)
    k = rand(UInt32, R)
    @simd for i in 0:N-1
        varphilox(i % UInt64, k, bits, bits)
    end
end
