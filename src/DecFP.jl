module DecFP
export Dec32, Dec64, Dec128, @d_str, @d32_str, @d64_str, @d128_str

using Compat

const libbid = joinpath(dirname(@__FILE__), "..", "deps", "libbid$WORD_SIZE")

const _buffer = Array(UInt8, 1024)

import Base.promote_rule
import Core.Intrinsics: box, unbox, bswap_int

# global pointers and dicts must be initialized at runtime (via __init__)
function __init__()
    global const rounding = cglobal((:__bid_IDEC_glbround, libbid), Cuint) # rounding mode
    global const flags = cglobal((:__bid_IDEC_glbflags, libbid), Cuint) # exception status
    unsafe_store!(flags, 0)

    # rounding modes, from bid_functions.h
    global const rounding_c2j = [RoundNearest, RoundDown, RoundUp, RoundToZero, RoundFromZero]
    global const rounding_j2c = [ rounding_c2j[i]=>Cuint(i-1) for i in 1:length(rounding_c2j) ]
end

# status flags from bid_functions.h:
const INVALID    = 0x01
const UNNORMAL   = 0x02
const DIVBYZERO  = 0x04
const OVERFLOW   = 0x08
const UNDERFLOW  = 0x10
const INEXACT    = 0x20

bidsym(w,s...) = string("__bid", w, "_", s...)

abstract DecimalFloatingPoint <: AbstractFloat
Base.get_rounding{T<:DecimalFloatingPoint}(::Type{T}) = rounding_c2j[unsafe_load(rounding)+1]
Base.set_rounding{T<:DecimalFloatingPoint}(::Type{T}, r::RoundingMode) = unsafe_store!(rounding, rounding_j2c[r])

for w in (32,64,128)
    BID = symbol(string("Dec",w))
    @eval bitstype $w $BID <: DecimalFloatingPoint
end

# quickly check whether s begins with "±nan"
function isnanstr(s::AbstractString)
    i = start(s)
    while !done(s, i)
        c, i = next(s, i)
        isspace(c) || break
    end
    done(s, i) && return false
    if (c == '+' || c == '-')
        c, i = next(s, i)
        done(s, i) && return false
    end
    lowercase(c) == 'n' || return false
    c, i = next(s, i)
    (!done(s, i) && lowercase(c) == 'a') || return false
    c, i = next(s, i)
    (done(s, i) && lowercase(c) == 'n') || return false
    return true
end

for w in (32,64,128)
    BID = symbol(string("Dec",w))
    T = eval(BID)
    Ti = eval(symbol(string("UInt",w)))

    # hack: we need an internal parsing function that doesn't check exceptions, since
    # flags isn't defined until __init__ runs.  Similarly for nextfloat/prevfloat
    @eval begin
        _parse(::Type{$BID}, s::AbstractString) =
            ccall(($(bidsym(w,"from_string")), libbid), $BID, (Ptr{UInt8},), s)
        _nextfloat(x::$BID) = ccall(($(bidsym(w,"nexttoward")), libbid), $BID, ($BID,Dec128), x, pinf128)
        _prevfloat(x::$BID) = ccall(($(bidsym(w,"nexttoward")), libbid), $BID, ($BID,Dec128), x, minf128)
        _sub(x::$BID, y::$BID) = ccall(($(bidsym(w,"sub")), libbid), $BID, ($BID,$BID), x, y)
    end

    @eval begin
        function Base.parse(::Type{$BID}, s::AbstractString)
            x = _parse($BID, s)
            if isnan(x) && !isnanstr(s)
                throw(ArgumentError("invalid number format $s"))
            end
            return xchk(x, InexactError)
        end

        function Base.show(io::IO, x::$BID)
            ccall(($(bidsym(w,"to_string")), libbid), Void, (Ptr{UInt8}, $BID), _buffer, x)
            write(io, pointer(_buffer), ccall(:strlen, Csize_t, (Ptr{UInt8},), _buffer))
        end

        Base.fma(x::$BID, y::$BID, z::$BID) = nox(ccall(($(bidsym(w,"fma")), libbid), $BID, ($BID,$BID,$BID), x, y, z))
        Base.muladd(x::$BID, y::$BID, z::$BID) = fma(x,y,z) # faster than x+y*z

        Base.one(::Union{Type{$BID},$BID}) = $(_parse(T, "1"))
        Base.zero(::Union{Type{$BID},$BID}) = $(_parse(T, "0"))

        Base.signbit(x::$BID) = $(zero(Ti)) != $(Ti(1) << (Ti(w - 1))) & reinterpret($Ti, x)
        Base.sign(x::$BID) = ifelse(signbit(x), $(_parse(T, "-1")), $(_parse(T, "1")))

        Base.nextfloat(x::$BID) = nox(_nextfloat(x))
        Base.prevfloat(x::$BID) = nox(_prevfloat(x))
        Base.eps(x::$BID) = ifelse(isfinite(x), xchk(nextfloat(x) - x, OVERFLOW), $(_parse(T, "NaN")))

        # the meaning of the exponent is different than for binary FP: it is 10^n, not 2^n:
        # Base.exponent(x::$BID) = nox(ccall(($(bidsym(w,"ilogb")), libbid), Cint, ($BID,), x))
        # Base.ldexp(x::$BID, n::Integer) = nox(ccall(($(bidsym(w,"ldexp")), libbid), $BID, ($BID,Cint), x, n))
    end

    for (f,c) in ((:isnan,"isNaN"), (:isinf,"isInf"), (:isfinite,"isFinite"), (:issubnormal,"isSubnormal"))
        @eval Base.$f(x::$BID) = ccall(($(bidsym(w,c)), libbid), Cint, ($BID,), x) != 0
    end

    for (f,c) in ((:+,"add"), (:-,"sub"), (:*,"mul"), (:/, "div"), (:hypot,"hypot"), (:atan2,"atan2"), (:^,"pow"), (:copysign,"copySign"))
        @eval Base.$f(x::$BID, y::$BID) = nox(ccall(($(bidsym(w,c)), libbid), $BID, ($BID,$BID), x, y))
    end

    for f in (:exp,:log,:sin,:cos,:tan,:asin,:acos,:atan,:sinh,:cosh,:tanh,:asinh,:acosh,:atanh,:log1p,:expm1,:log10,:log2,:exp2,:exp10,:erf,:erfc,:lgamma,:sqrt,:cbrt,:abs)
        @eval Base.$f(x::$BID) = xchk(ccall(($(bidsym(w,f)), libbid), $BID, ($BID,), x), INVALID)
    end
    for (f,c) in ((:gamma,"tgamma"), (:-,"negate"), (:round,"nearbyint"))
        @eval Base.$f(x::$BID) = xchk(ccall(($(bidsym(w,c)), libbid), $BID, ($BID,), x), INVALID)
    end

    for (f,c) in ((:(==),"quiet_equal"), (:>,"quiet_greater"), (:<,"quiet_less"), (:(>=), "quiet_greater_equal"), (:(<=), "quiet_less_equal"))
        @eval Base.$f(x::$BID, y::$BID) = nox(ccall(($(bidsym(w,c)), libbid), Cint, ($BID,$BID), x, y) != 0)
    end

    for Tf in (Float32,Float64)
        bT = string("binary",sizeof(Tf)*8)
        @eval begin
            Base.convert(::Type{$Tf}, x::$BID) = nox(ccall(($(bidsym(w,"to_",bT)), libbid), $Tf, ($BID,), x))
            Base.convert(::Type{$BID}, x::$Tf) = nox(ccall(($(string("__",bT,"_to_","bid",w)), libbid), $BID, ($Tf,), x))
        end
    end

    for c in (:π, :e, :γ, :catalan, :φ)
        @eval begin
            Base.convert(::Type{$BID}, ::Irrational{$(QuoteNode(c))}) = $(_parse(T, with_bigfloat_precision(256) do
                                                                                      string(BigFloat(eval(c)))
                                                                                  end))
            promote_rule(::Type{$BID}, ::Type{Irrational}) = $BID
        end
    end

    for w′ in (32,64,128)
        BID′ = symbol(string("Dec",w′))
        if w > w′
            @eval promote_rule(::Type{$BID}, ::Type{$BID′}) = $BID
        end
        if w != w′
            @eval Base.convert(::Type{$BID}, x::$BID′) = xchk(ccall(($(string("__bid",w′,"_to_","bid",w)), libbid), $BID, ($BID′,), x), INEXACT)
        end

        # promote binary*decimal -> decimal, for consistency with other operations above
        # (there doesn't seem to be any clear standard for this)
        if w′ <= 64
            FP′ = symbol(string("Float",w′))
            @eval promote_rule(::Type{$BID}, ::Type{$FP′}) = $(symbol(string("Dec",max(w,w′))))
            for i′ in ("Int$w′", "UInt$w′")
                Ti′ = eval(symbol(i′))
                @eval begin
                    Base.convert(::Type{$BID}, x::$Ti′) = nox(ccall(($(bidsym(w,"from_",lowercase(i′))), libbid), $BID, ($Ti′,), x))
                end
            end
        end
    end

    for w′ in (8,16,32,64)
        for i′ in ("Int$w′", "UInt$w′")
            Ti′ = eval(symbol(i′))
            @eval begin
                Base.floor(::Type{$Ti′}, x::$BID) = xchk(ccall(($(bidsym(w,"to_",lowercase(i′),"_xfloor")), libbid), $Ti′, ($BID,), x), InexactError, INVALID | OVERFLOW)
                Base.ceil(::Type{$Ti′}, x::$BID) = xchk(ccall(($(bidsym(w,"to_",lowercase(i′),"_xceil")), libbid), $Ti′, ($BID,), x), InexactError, INVALID | OVERFLOW)
                Base.convert(::Type{$Ti′}, x::$BID) = xchk(ccall(($(bidsym(w,"to_",lowercase(i′),"_xfloor")), libbid), $Ti′, ($BID,), x), InexactError)
            end
        end
    end
    
    @eval Base.bswap(x::$BID) = reinterpret($BID, bswap(reinterpret($Ti, x)))
    @eval Base.convert(::Type{Float16}, x::$BID) = convert(Float16, convert(Float32, x))
end # widths w

# the complex-sqrt function in base doesn't work for use, because it requires base-2 ldexp
function Base.sqrt{T<:DecimalFloatingPoint}(z::Complex{T})
    x, y = reim(z)
    x==y==0 && return Complex(zero(x),y)
    ρ = sqrt((abs(x) + hypot(x,y)) * 0.5)
    ξ = ρ
    η = y
    if isfinite(η) η=(η/ρ)/2 end
    if x<0
        ξ = abs(η)
        η = copysign(ρ,y)
    end
    Complex(ξ,η)
end

# used for next/prevfloat:
const pinf128 = _parse(Dec128, "+Inf")
const minf128 = _parse(Dec128, "-Inf")

for T in (Dec32,Dec64,Dec128)
    @eval begin
        Base.eps(::Type{$T}) = $(_sub(_nextfloat(one(T)), one(T)))
        Base.typemax(::Type{$T}) = $(_parse(T, "+inf"))
        Base.typemin(::Type{$T}) = $(_parse(T, "-inf"))
        Base.realmax(::Type{$T}) = $(_prevfloat(_parse(T, "+inf")))
        Base.realmin(::Type{$T}) = $(_nextfloat(zero(T)))
    end
end

Base.convert{F<:DecimalFloatingPoint}(T::Type{F}, x::Union{Int8,UInt8,Int16,UInt16}) = F(Int32(x))
Base.convert{F<:DecimalFloatingPoint}(T::Type{F}, x::Float16) = F(Float32(x))
promote_rule{F<:DecimalFloatingPoint}(::Type{F}, ::Type{Float16}) = F
promote_rule{F<:DecimalFloatingPoint,T<:Union{Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64}}(::Type{F}, ::Type{T}) = F

# so that mathconsts get promoted to Dec32, not Dec64, like Float32
promote_rule{s,F<:DecimalFloatingPoint}(::Type{Irrational{s}}, ::Type{F}) = F
promote_rule{s,F<:DecimalFloatingPoint}(::Type{Irrational{s}}, T::Type{Complex{F}}) = T

macro d_str(s, flags...) parse(Dec64, s) end
macro d32_str(s, flags...) parse(Dec32, s) end
macro d64_str(s, flags...) parse(Dec64, s) end
macro d128_str(s, flags...) parse(Dec128, s) end

# clear exception flags and return x
function nox(x)
    unsafe_store!(flags, 0)
    return x
end

# check exception flags in mask & throw, otherwise returning x;
# always clearing exceptions
function xchk(x, mask::Integer=0x3f)
    f = unsafe_load(flags)
    unsafe_store!(flags, 0)
    if f & mask != 0
        f & INEXACT != 0 && throw(InexactError())
        f & OVERFLOW != 0 && throw(OverflowError())
        f & DIVBYZERO != 0 && throw(DivideError())
        f & INVALID != 0 && throw(DomainError())
        f & UNDERFLOW != 0 && error("underflow")
        f & UNNORMAL != 0 && error("unnormal")
    end
    return x
end

function xchk{E<:Exception}(x, exc::Type{E}, mask::Integer=0x3f)
    f = unsafe_load(flags)
    unsafe_store!(flags, 0)
    f & mask != 0 && throw(exc())
    return x
end

end # module
