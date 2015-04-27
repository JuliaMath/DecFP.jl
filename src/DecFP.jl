module DecFP
export Dec32, Dec64, Dec128, @d_str, @d32_str, @d64_str, @d128_str

using Compat

const libbid = joinpath(dirname(@__FILE__), "..", "deps", "libbid$WORD_SIZE")

const _buffer = Array(UInt8, 1024)

import Base.promote_rule

# global pointers and dicts must be initialized at runtime (via __init__)
function __init__()
    global const rounding = cglobal((:__bid_IDEC_glbround, libbid), Cuint) # rounding mode
    global const flags = cglobal((:__bid_IDEC_glbflags, libbid), Cuint) # exception status

    # rounding modes, from bid_functions.h
    global const rounding_c2j = [RoundNearest, RoundDown, RoundUp, RoundToZero, RoundFromZero]
    global const rounding_j2c = [ rounding_c2j[i]=>Cuint(i-1) for i in 1:length(rounding_c2j) ]
end

bidsym(w,s...) = string("__bid", w, "_", s...)

for w in (32,64,128)
    BID = symbol(string("Dec",w))
    @eval begin
        bitstype $w $BID <: FloatingPoint

        function Base.parse(::Type{$BID}, s::AbstractString)
            x = ccall(($(bidsym(w,"from_string")), libbid), $BID, (Ptr{UInt8},), s)
            if isnan(x)
                # fixme: check whether s is "nan" etc.
                throw(ArgumentError("invalid number format $s"))
            end
            return x
        end

        function Base.show(io::IO, x::$BID)
            ccall(($(bidsym(w,"to_string")), libbid), Void, (Ptr{UInt8}, $BID), _buffer, x)
            write(io, pointer(_buffer), ccall(:strlen, Csize_t, (Ptr{UInt8},), _buffer))
        end

        Base.get_rounding(::Type{$BID}) = rounding_c2j[unsafe_load(rounding)+1]
        Base.set_rounding(::Type{$BID}, r::RoundingMode) = unsafe_store!(rounding, rounding_j2c[r])

        Base.fma(x::$BID, y::$BID, z::$BID) = ccall(($(bidsym(w,"fma")), libbid), $BID, ($BID,$BID,$BID), x, y, z)
        Base.muladd(x::$BID, y::$BID, z::$BID) = fma(x,y,z) # faster than x+y*z
        Base.gamma(x::$BID) = ccall(($(bidsym(w,"tgamma")), libbid), $BID, ($BID,), x)
        Base.$(:-)(x::$BID) = ccall(($(bidsym(w,"negate")), libbid), $BID, ($BID,), x)
    end

    for (f,c) in ((:isnan,"isNaN"), (:isinf,"isInf"), (:isfinite,"isFinite"), (:issubnormal,"isSubnormal"))
        @eval Base.$f(x::$BID) = ccall(($(bidsym(w,c)), libbid), Cint, ($BID,), x) != 0
    end

    for (f,c) in ((:+,"add"), (:-,"sub"), (:*,"mul"), (:/, "div"), (:hypot,"hypot"), (:atan2,"atan2"), (:mod,"fmod"), (:^,"pow"), (:copysign,"copySign"))
        @eval Base.$f(x::$BID, y::$BID) = ccall(($(bidsym(w,c)), libbid), $BID, ($BID,$BID), x, y)
    end

    for f in (:exp,:log,:sin,:cos,:tan,:asin,:acos,:atan,:sinh,:cosh,:tanh,:asinh,:acosh,:atanh,:log1p,:expm1,:log10,:log2,:exp2,:exp10,:erf,:erfc,:lgamma,:sqrt,:cbrt,:abs)
        @eval Base.$f(x::$BID) = ccall(($(bidsym(w,f)), libbid), $BID, ($BID,), x)
    end

    for (f,c) in ((:(==),"quiet_equal"), (:>,"quiet_greater"), (:<,"quiet_less"), (:(>=), "quiet_greater_equal"), (:(<=), "quiet_less_equal"))
        @eval Base.$f(x::$BID, y::$BID) = ccall(($(bidsym(w,c)), libbid), Cint, ($BID,$BID), x, y) != 0
    end

    for T in (Float32,Float64)
        bT = string("binary",sizeof(T)*8)
        @eval begin
            Base.convert(::Type{$T}, x::$BID) = ccall(($(bidsym(w,"to_",bT)), libbid), $T, ($BID,), x)
            Base.convert(::Type{$BID}, x::$T) = ccall(($(string("__",bT,"_to_","bid",w)), libbid), $BID, ($T,), x)
        end
    end

    @eval begin
        const $(symbol(string("one",w))) = parse($BID, "1")
        const $(symbol(string("zero",w))) = parse($BID, "0")
        Base.one(::Union(Type{$BID},$BID)) = $(symbol(string("one",w)))
        Base.zero(::Union(Type{$BID},$BID)) = $(symbol(string("zero",w)))
    end

    for c in (:π, :e, :γ, :catalan, :φ)
        @eval begin
            const $(symbol(string(c,w))) = parse($BID, with_bigfloat_precision(256) do
                                                           string(BigFloat($c))
                                                       end)
            Base.convert(::Type{$BID}, ::MathConst{$(QuoteNode(c))}) = $(symbol(string(c,w)))
            promote_rule(::Type{$BID}, ::Type{MathConst}) = $BID
        end
    end

    for w′ in (32,64,128)
        BID′ = symbol(string("Dec",w′))
        if w > w′
            @eval promote_rule(::Type{$BID}, ::Type{$BID′}) = $BID
        end

        # promote binary*decimal -> decimal, for consistency with other operations above
        # (there doesn't seem to be any clear standard for this)
        if w′ <= 64
            FP′ = symbol(string("Float",w′))
            @eval promote_rule(::Type{$BID}, ::Type{$FP′}) = $(symbol(string("Dec",max(w,w′))))
        end
    end

    @eval promote_rule{T<:Union(Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64,Int128,UInt128)}(::Type{$BID}, ::Type{T}) = $BID
end # widths w

macro d_str(s, flags...) parse(Dec64, s) end
macro d32_str(s, flags...) parse(Dec32, s) end
macro d64_str(s, flags...) parse(Dec64, s) end
macro d128_str(s, flags...) parse(Dec128, s) end

end # module
