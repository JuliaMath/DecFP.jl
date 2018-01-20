__precompile__(true)
module DecFP

using Compat, Compat.Printf, Compat.Unicode

export Dec32, Dec64, Dec128, @d_str, @d32_str, @d64_str, @d128_str

const libbid = joinpath(dirname(@__FILE__), "..", "deps", "libbid$(Sys.WORD_SIZE)")

const _buffer = fill(0x00, 1024)

import Base.promote_rule
import Base.Grisu.DIGITS

const rounding = Ref{Ptr{Cuint}}()
const flags = Ref{Ptr{Cuint}}()
# rounding modes, from bid_functions.h

const rounding_c2j = [RoundNearest, RoundDown, RoundUp, RoundToZero, RoundFromZero]
const rounding_j2c = Dict{RoundingMode, UInt32}([(rounding_c2j[i], Cuint(i-1)) for i in 1:length(rounding_c2j)])

# global pointers and dicts must be initialized at runtime (via __init__)
function __init__()
    global rounding[] = cglobal((:__bid_IDEC_glbround, libbid), Cuint) # rounding mode
    global flags[] = cglobal((:__bid_IDEC_glbflags, libbid), Cuint) # exception status
    unsafe_store!(flags[], 0)
end

# status flags from bid_functions.h:
const INVALID    = 0x01
const UNNORMAL   = 0x02
const DIVBYZERO  = 0x04
const OVERFLOW   = 0x08
const UNDERFLOW  = 0x10
const INEXACT    = 0x20

bidsym(w,s...) = string("__bid", w, "_", s...)

abstract type DecimalFloatingPoint <: AbstractFloat end
Base.rounding(::Type{T}) where {T<:DecimalFloatingPoint} = rounding_c2j[unsafe_load(rounding[])+1]
Base.setrounding(::Type{T}, r::RoundingMode) where {T<:DecimalFloatingPoint} = unsafe_store!(rounding[], rounding_j2c[r])

for w in (32,64,128)
    BID = Symbol(string("Dec",w))
    Ti = Symbol(string("UInt",w))
    @eval struct $BID <: DecimalFloatingPoint
        x::$Ti
        $BID(x) = convert($BID, x)
        Base.reinterpret(::Type{$BID}, x::$Ti) = new(x)
    end
    @eval $BID(x::Rational{T}) where {T} = convert($BID, x)
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
    (c == 'n' || c == 'N') || return false
    c, i = next(s, i)
    (!done(s, i) && (c == 'a' || c == 'A')) || return false
    c, i = next(s, i)
    (done(s, i) && (c == 'n' || c == 'N')) || return false
    return true
end

function Base.show(io::IO, x::DecimalFloatingPoint)
    s = @sprintf("%g", x)
    if ismatch(r"^-?\d+$", s)
        s *= ".0"
    end
    print(io, s)
end

for w in (32,64,128)
    BID = Symbol(string("Dec",w))
    Ti = eval(Symbol(string("UInt",w)))
    T = eval(BID)

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
            return xchk(x, InexactError, :parse, $BID, s)
        end

        $BID(x::AbstractString) = parse($BID, x)

        function Base.Printf.fix_dec(x::$BID, n::Int)
            if n > length(DIGITS) - 1
                n = length(DIGITS) - 1
            end
            # rounded = round(x * exp10($BID(n)), RoundNearestTiesAway)
            rounded = xchk(ccall(($(bidsym(w,"round_integral_nearest_away")), libbid), $BID, ($BID,), x * exp10($BID(n))), InexactError, :round, $BID, x, mask=INVALID | OVERFLOW)
            if rounded == 0
                DIGITS[1] = UInt8('0')
                return Int32(1), Int32(1), signbit(x)
            end
            ccall(($(bidsym(w,"to_string")), libbid), Cvoid, (Ptr{UInt8}, $BID), _buffer, rounded)
            trailing_zeros = 0
            i = 2
            while _buffer[i] != UInt8('E')
                DIGITS[i - 1] = _buffer[i]
                if _buffer[i] == UInt8('0')
                    trailing_zeros += 1
                else
                    trailing_zeros = 0
                end
                i += 1
            end
            ndigits = i - 2
            len = ndigits - trailing_zeros
            i += 1
            if _buffer[i] == UInt8('+')
                expsign = +1
            elseif _buffer[i] == UInt8('-')
                expsign = -1
            end
            exponent = 0
            i += 1
            while _buffer[i] != 0x00
                exponent = exponent * 10 + _buffer[i] - UInt8('0')
                i += 1
            end
            exponent *= expsign
            pt = ndigits + exponent - n
            neg = signbit(x)
            return Int32(len), Int32(pt), neg
        end

        function Base.Printf.ini_dec(x::$BID, n::Int)
            if n > length(DIGITS) - 1
                n = length(DIGITS) - 1
            end
            if x == 0
                for i = 1:n
                    DIGITS[i] = UInt8('0')
                end
                return Int32(1), Int32(1), signbit(x)
            end
            normalized_exponent = nox(ccall(($(bidsym(w,"ilogb")), libbid), Cint, ($BID,), x))
            # rounded = round(x * exp10($BID(n - 1 - normalized_exponent)), RoundNearestTiesAway)
            rounded = xchk(ccall(($(bidsym(w,"round_integral_nearest_away")), libbid), $BID, ($BID,), x * exp10($BID(n - 1 - normalized_exponent))), InexactError, :round, $BID, x, mask=INVALID | OVERFLOW)
            rounded_exponent = nox(ccall(($(bidsym(w,"ilogb")), libbid), Cint, ($BID,), rounded))
            ccall(($(bidsym(w,"to_string")), libbid), Cvoid, (Ptr{UInt8}, $BID), _buffer, rounded)
            i = 2
            while _buffer[i] != UInt8('E')
                DIGITS[i - 1] = _buffer[i]
                i += 1
            end
            pt = normalized_exponent + rounded_exponent - n + 2
            neg = signbit(x)
            return Int32(n), Int32(pt), neg
        end

        Base.fma(x::$BID, y::$BID, z::$BID) = nox(ccall(($(bidsym(w,"fma")), libbid), $BID, ($BID,$BID,$BID), x, y, z))
        Base.muladd(x::$BID, y::$BID, z::$BID) = fma(x,y,z) # faster than x+y*z

        Base.one(::Union{Type{$BID},$BID}) = $(_parse(T, "1"))
        Base.zero(::Union{Type{$BID},$BID}) = $(_parse(T, "0"))

        Base.signbit(x::$BID) = $(zero(Ti)) != $(Ti(1) << (Ti(w - 1))) & x.x
        Base.sign(x::$BID) = ifelse(signbit(x), $(_parse(T, "-1")), $(_parse(T, "1")))

        Base.nextfloat(x::$BID) = nox(_nextfloat(x))
        Base.prevfloat(x::$BID) = nox(_prevfloat(x))
        Base.eps(x::$BID) = ifelse(isfinite(x), xchk(nextfloat(x) - x, OVERFLOW, "$($BID) value overflow"), $(_parse(T, "NaN")))

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

    for f in (:exp,:log,:sin,:cos,:tan,:asin,:acos,:atan,:sinh,:cosh,:tanh,:asinh,:acosh,:atanh,:log1p,:expm1,:log10,:log2,:exp2,:exp10,:lgamma,:sqrt,:cbrt,:abs)
        @eval Base.$f(x::$BID) = xchk(ccall(($(bidsym(w,f)), libbid), $BID, ($BID,), x), "invalid operation '$($f)' on $($BID)", mask=INVALID)
    end

    for (f,c) in ((:gamma,"tgamma"), (:-,"negate"), (:round,"nearbyint"))
        @eval Base.$f(x::$BID) = xchk(ccall(($(bidsym(w,c)), libbid), $BID, ($BID,), x), "invalid operation '$($c)' on $($BID)", mask=INVALID)
    end

    for (f,c) in ((:(==),"quiet_equal"), (:>,"quiet_greater"), (:<,"quiet_less"), (:(>=), "quiet_greater_equal"), (:(<=), "quiet_less_equal"))
        @eval Base.$f(x::$BID, y::$BID) = nox(ccall(($(bidsym(w,c)), libbid), Cint, ($BID,$BID), x, y) != 0)
    end

    for Tf in (Float32,Float64)
        bT = string("binary",sizeof(Tf)*8)
        @eval begin
            Base.convert(::Type{$Tf}, x::$BID) = nox(ccall(($(bidsym(w,"to_",bT)), libbid), $Tf, ($BID,), x))
            Base.$(Symbol("$Tf"))(x::$BID) = convert($Tf, x)
            Base.convert(::Type{$BID}, x::$Tf) = nox(ccall(($(string("__",bT,"_to_","bid",w)), libbid), $BID, ($Tf,), x))
        end
    end

    for c in (:π, :e, :ℯ, :γ, :catalan, :φ)
        @eval begin
            Base.convert(::Type{$BID}, ::Irrational{$(QuoteNode(c))}) = $(_parse(T, setprecision(256) do
                                                                                      string(BigFloat(isdefined(Base, :MathConstants) ? eval(Base.MathConstants, c) : eval(c)))
                                                                                  end))
        end
    end

    @eval promote_rule(::Type{$BID}, ::Type{Irrational}) = $BID

    for w′ in (32,64,128)
        BID′ = Symbol(string("Dec",w′))
        if w > w′
            @eval promote_rule(::Type{$BID}, ::Type{$BID′}) = $BID
        end
        if w != w′
            @eval Base.convert(::Type{$BID}, x::$BID′) = xchk(ccall(($(string("__bid",w′,"_to_","bid",w)), libbid), $BID, ($BID′,), x), INEXACT, :convert, $BID, x)
        end

        # promote binary*decimal -> decimal, for consistency with other operations above
        # (there doesn't seem to be any clear standard for this)
        if w′ <= 64
            FP′ = Symbol(string("Float",w′))
            @eval promote_rule(::Type{$BID}, ::Type{$FP′}) = $(Symbol(string("Dec",max(w,w′))))
            for (i′, i′str) in (("Int$w′", "int$w′"), ("UInt$w′", "uint$w′"))
                Ti′ = eval(Symbol(i′))
                @eval begin
                    Base.convert(::Type{$BID}, x::$Ti′) = nox(ccall(($(bidsym(w,"from_",i′str)), libbid), $BID, ($Ti′,), x))
                end
            end
        end
    end

    for w′ in (8,16,32,64)
        for (i′, i′str) in (("Int$w′", "int$w′"), ("UInt$w′", "uint$w′"))
            Ti′ = eval(Symbol(i′))
            @eval begin
                Base.trunc(::Type{$Ti′}, x::$BID) = xchk(ccall(($(bidsym(w,"to_",i′str,"_xint")), libbid), $Ti′, ($BID,), x), InexactError, :trunc, $BID, x, mask=INVALID | OVERFLOW)
                Base.floor(::Type{$Ti′}, x::$BID) = xchk(ccall(($(bidsym(w,"to_",i′str,"_xfloor")), libbid), $Ti′, ($BID,), x), InexactError, :floor, $BID, x, mask=INVALID | OVERFLOW)
                Base.ceil(::Type{$Ti′}, x::$BID) = xchk(ccall(($(bidsym(w,"to_",i′str,"_xceil")), libbid), $Ti′, ($BID,), x), InexactError, :ceil, $BID, x, mask=INVALID | OVERFLOW)
                Base.round(::Type{$Ti′}, x::$BID) = xchk(ccall(($(bidsym(w,"to_",i′str,"_xrnint")), libbid), $Ti′, ($BID,), x), InexactError, :round, $BID, x, mask=INVALID | OVERFLOW)
                Base.round(::Type{$Ti′}, x::$BID, ::RoundingMode{:NearestTiesAway}) = xchk(ccall(($(bidsym(w,"to_",i′str,"_xrninta")), libbid), $Ti′, ($BID,), x), InexactError, :round, $BID, x, mask=INVALID | OVERFLOW)
                Base.convert(::Type{$Ti′}, x::$BID) = xchk(ccall(($(bidsym(w,"to_",i′str,"_xfloor")), libbid), $Ti′, ($BID,), x), InexactError, :convert, $BID, x)
                Base.$(Symbol("$Ti′"))(x::$BID) = convert($Ti′, x)
            end
        end
    end

    @eval Base.bswap(x::$BID) = reinterpret($BID, bswap(x.x))
    @eval Base.convert(::Type{Float16}, x::$BID) = convert(Float16, convert(Float32, x))
    @eval Base.Float16(x::$BID) = convert(Float16, x)
    @eval Base.reinterpret(::Type{$Ti}, x::$BID) = x.x
end # widths w

Base.trunc(::Type{Integer}, x::DecimalFloatingPoint) = trunc(Int, x)
Base.floor(::Type{Integer}, x::DecimalFloatingPoint) = floor(Int, x)
Base.ceil(::Type{Integer}, x::DecimalFloatingPoint) = ceil(Int, x)
Base.round(::Type{Integer}, x::DecimalFloatingPoint) = round(Int, x)
Base.round(::Type{Integer}, x::DecimalFloatingPoint, ::RoundingMode{:NearestTiesAway}) = round(Int, x, RoundNearestTiesAway)
Base.convert(::Type{Integer}, x::DecimalFloatingPoint) = convert(Int, x)

Base.round(::Type{T}, x::DecimalFloatingPoint, ::RoundingMode{:Nearest}) where {T<:Integer} = round(T, x)
function Base.round(::Type{T}, x::DecimalFloatingPoint, ::RoundingMode{:NearestTiesUp}) where {T<:Integer}
    y = floor(T, x)
    ifelse(x==y, y, copysign(floor(T, 2*x-y), x))
end
Base.round(::Type{T}, x::DecimalFloatingPoint, ::RoundingMode{:ToZero}) where {T<:Integer} = trunc(T, x)
Base.round(::Type{T}, x::DecimalFloatingPoint, ::RoundingMode{:FromZero}) where {T<:Integer} = (x>=0 ? ceil(T, x) : floor(T, x))
Base.round(::Type{T}, x::DecimalFloatingPoint, ::RoundingMode{:Up}) where {T<:Integer} = ceil(T, x)
Base.round(::Type{T}, x::DecimalFloatingPoint, ::RoundingMode{:Down}) where {T<:Integer} = floor(T, x)

# the complex-sqrt function in base doesn't work for use, because it requires base-2 ldexp
function Base.sqrt(z::Complex{T}) where {T<:DecimalFloatingPoint}
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

Base.convert(T::Type{F}, x::Union{Int8,UInt8,Int16,UInt16}) where {F<:DecimalFloatingPoint} = F(Int32(x))
Base.convert(T::Type{F}, x::Integer) where {F<:DecimalFloatingPoint} = F(Int64(x))
Base.convert(T::Type{F}, x::Unsigned) where {F<:DecimalFloatingPoint} = F(UInt64(x))
Base.convert(T::Type{F}, x::Rational) where {F<:DecimalFloatingPoint} = F(x.num) / F(x.den)
Base.convert(T::Type{F}, x::Float16) where {F<:DecimalFloatingPoint} = F(Float32(x))
promote_rule(::Type{F}, ::Type{Float16}) where {F<:DecimalFloatingPoint} = F
promote_rule(::Type{F}, ::Type{T}) where {F<:DecimalFloatingPoint,T<:Union{Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64}} = F

# so that mathconsts get promoted to Dec32, not Dec64, like Float32
promote_rule(::Type{Irrational{s}}, ::Type{F}) where {s,F<:DecimalFloatingPoint} = F
promote_rule(::Type{Irrational{s}}, T::Type{Complex{F}}) where {s,F<:DecimalFloatingPoint} = T

macro d_str(s, flags...) parse(Dec64, s) end
macro d32_str(s, flags...) parse(Dec32, s) end
macro d64_str(s, flags...) parse(Dec64, s) end
macro d128_str(s, flags...) parse(Dec128, s) end

# clear exception flags and return x
function nox(x)
    unsafe_store!(flags[], 0)
    return x
end

# check exception flags in mask & throw, otherwise returning x;
# always clearing exceptions
function xchk(x, args...; mask::Integer=0x3f)
    f = unsafe_load(flags[])
    unsafe_store!(flags[], 0)
    if f & mask != 0
        f & INEXACT != 0 && throw(InexactError(args...))
        f & OVERFLOW != 0 && throw(OverflowError(args...))
        f & DIVBYZERO != 0 && throw(DivideError())
        f & INVALID != 0 && throw(DomainError(args...))
        f & UNDERFLOW != 0 && error("underflow")
        f & UNNORMAL != 0 && error("unnormal")
    end
    return x
end

function xchk(x, exc::Type{E}, args...; mask::Integer=0x3f) where {E<:Exception}
    f = unsafe_load(flags[])
    unsafe_store!(flags[], 0)
    f & mask != 0 && throw(exc(args...))
    return x
end

end # module
