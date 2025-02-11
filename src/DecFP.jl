
module DecFP

using DecFP_jll

import Printf, Random, SpecialFunctions

export Dec32, Dec64, Dec128, @d_str, @d32_str, @d64_str, @d128_str, exponent10, ldexp10, sigexp

import Base.promote_rule

# Check exception flags in mask & throw, otherwise returning x;
# This is a macros so that the error message is only evaluated if an exception occurs.
macro xchk(x, flags, exc, args...)
    mask=0x3f
    if !isempty(args) && Meta.isexpr(args[end], :(=)) && args[end].args[1] == :mask # mask=... keyword at end
        mask = esc(args[end].args[2])
        args = args[1:end-1]
    end
    quote
        ret = $(esc(x))
        if !isnothing($exc)
            $(esc(flags)) & $mask != 0 && throw($exc($(map(esc,args)...)))
        end
        ret
    end
end

#############################################################################

@enum DecFPRoundingMode begin
    DecFPRoundNearest
    DecFPRoundDown
    DecFPRoundUp
    DecFPRoundToZero
    DecFPRoundFromZero
end

Base.convert(::Type{DecFPRoundingMode}, ::RoundingMode{:Nearest})  = DecFPRoundNearest
Base.convert(::Type{DecFPRoundingMode}, ::RoundingMode{:Down})     = DecFPRoundDown
Base.convert(::Type{DecFPRoundingMode}, ::RoundingMode{:Up})       = DecFPRoundUp
Base.convert(::Type{DecFPRoundingMode}, ::RoundingMode{:ToZero})   = DecFPRoundToZero
Base.convert(::Type{DecFPRoundingMode}, ::RoundingMode{:FromZero}) = DecFPRoundFromZero

function Base.convert(::Type{RoundingMode}, r::DecFPRoundingMode)
    if r == DecFPRoundNearest
        return RoundNearest
    elseif r == DecFPRoundDown
        return RoundDown
    elseif r == DecFPRoundUp
        return RoundUp
    elseif r == DecFPRoundToZero
        return RoundToZero
    elseif r == DecFPRoundFromZero
        return RoundFromZero
    else
        throw(ArgumentError("invalid DecFP rounding mode code: $c"))
    end
end

# internal task-local buffer for I/O and string conversions
const _STRINGBUFFER_KEY = :DecFP_stringbuffer_abb78e082af23329 # unique key
_stringbuffer() = get!(() -> zeros(UInt8, 1024), task_local_storage(), _STRINGBUFFER_KEY)::Vector{UInt8}

# status flags from bid_functions.h:
const INVALID    = 0x01
const UNNORMAL   = 0x02
const DIVBYZERO  = 0x04
const OVERFLOW   = 0x08
const UNDERFLOW  = 0x10
const INEXACT    = 0x20

bidsym(w,s...) = string("__bid", w, "_", s...)

abstract type DecimalFloatingPoint <: AbstractFloat end

# fixme: make per-type? is task-local worth it or should we just make it global?
const _ROUNDINGMODE_KEY = :DecFP_roundingmode_abb78e082af23329 # unique key
_roundingmode() = get!(task_local_storage(), _ROUNDINGMODE_KEY, DecFPRoundNearest)::DecFPRoundingMode
Base.Rounding.rounding_raw(::Type{T}) where {T<:DecimalFloatingPoint} = _roundingmode()
Base.Rounding.setrounding_raw(::Type{T}, r::DecFPRoundingMode) where {T<:DecimalFloatingPoint} =
    task_local_storage(_ROUNDINGMODE_KEY, r)
Base.Rounding.setrounding_raw(f::Function, ::Type{T}, r::DecFPRoundingMode) where {T<:DecimalFloatingPoint} =
    task_local_storage(f, _ROUNDINGMODE_KEY, r)

Base.Rounding.rounding(::Type{T}) where {T<:DecimalFloatingPoint} =
    convert(RoundingMode, Base.Rounding.rounding_raw(T))
Base.Rounding.setrounding(::Type{T}, r::RoundingMode) where {T<:DecimalFloatingPoint} =
    Base.Rounding.setrounding_raw(T, convert(DecFPRoundingMode, r))
Base.Rounding.setrounding(f::Function, ::Type{T}, r::RoundingMode) where {T<:DecimalFloatingPoint} =
    Base.Rounding.setrounding_raw(f, T, convert(DecFPRoundingMode, r))

@static if Sys.ARCH == :arm64 || Sys.ARCH == :aarch64 || (isdefined(Base, :BinaryPlatforms) && Base.BinaryPlatforms.arch(Base.BinaryPlatforms.HostPlatform()) == "aarch64")
    # primitive types aren't working yet on ARM64 for some reason?
    struct Dec32 <: DecimalFloatingPoint
        x::UInt32
        Base.reinterpret(::Type{Dec32}, x::UInt32) = new(x)
    end
    Base.reinterpret(::Type{UInt32}, x::Dec32) = x.x
    Base.bswap(x::Dec32) = reinterpret(Dec32, bswap(x.x))
    struct Dec64 <: DecimalFloatingPoint
        x::UInt64
        Base.reinterpret(::Type{Dec64}, x::UInt64) = new(x)
    end
    Base.reinterpret(::Type{UInt64}, x::Dec64) = x.x
    Base.bswap(x::Dec64) = reinterpret(Dec64, bswap(x.x))
else
    primitive type Dec32 <: DecimalFloatingPoint 32 end
    primitive type Dec64 <: DecimalFloatingPoint 64 end
    Base.bswap(x::Dec32) = reinterpret(Dec32, bswap(reinterpret(UInt32, x)))
    Base.bswap(x::Dec64) = reinterpret(Dec64, bswap(reinterpret(UInt64, x)))
end
struct Dec128 <: DecimalFloatingPoint
    x::UInt128
    Base.reinterpret(::Type{Dec128}, x::UInt128) = new(x)
end
Base.reinterpret(::Type{UInt128}, x::Dec128) = x.x
Base.bswap(x::Dec128) = reinterpret(Dec128, bswap(x.x))

for w in (32,64,128)
    BID = Symbol(string("Dec",w))
    Ti = Symbol(string("UInt",w))

    @eval function $BID(x::Real, mode::RoundingMode)
        setrounding($BID, mode) do
            convert($BID, x)
        end
    end

    @eval function $BID(sign::Integer, significand::Integer, exponent::Integer)
        -1 <= sign <= 1 || throw(DomainError(sign, "sign must be -1, 0, or +1"))
        significand == 0 && return flipsign(zero($BID), sign)
        sign == 0 && throw(DomainError(sign, "sign must be -1 or +1 for non-zero significand"))
        p = 9 * $w ÷ 32 - 2
        emax = 3 * 2^($w ÷ 16 + 3)
        bias = emax + p - 2
        bemax = emax * 2 - 1
        ebits = $w ÷ 16 + 6
        sbits = $w - ebits - 1
        sb = signbit(sign) ? one($Ti) << ($w - 1) : zero($Ti)
        e = exponent + bias
        s = $Ti(abs(significand))
        while s >= _int_maxintfloat($BID)
            q, r = divrem(s, $Ti(10))
            r != 0 && throw(InexactError(Symbol($BID), $BID, (sign, significand, exponent)))
            s = q
            e += 1
        end
        while e > bemax && s < _int_maxintfloat($BID) ÷ $Ti(10)
            s *= $Ti(10)
            e -= 1
        end
        e > bemax && throw(InexactError(Symbol($BID), $BID, (sign, significand, exponent)))
        while e < 0
            q, r = divrem(s, $Ti(10))
            r != 0 && throw(InexactError(Symbol($BID), $BID, (sign, significand, exponent)))
            s = q
            e += 1
        end
        if s < one($Ti) << sbits
            return reinterpret($BID, sb | ($Ti(e) << sbits) | s)
        end
        return reinterpret($BID, sb | ($Ti(0x3) << ($w - 3)) | ($Ti(e) << (sbits - 2)) | (s & (typemax($Ti) >> (ebits + 3))))
    end

    @eval $BID(significand::Integer, exponent::Integer) = $BID(sign(significand), significand, exponent)

    @eval @doc """
        $($BID)(x::Union{Real, AbstractString} [, mode::RoundingMode])
        $($BID)([sign::Integer,] significand::Integer, exponent::Integer)

    Create a $($w)-bit IEEE 754-2008 decimal floating point number. The `mode` argument
    specifies the direction in which the result should be rounded if the conversion cannot
    be done exactly. If not provided, the `mode` is set by the current
    `rounding(DecFP.DecimalFloatingPoint)` mode.

    `$($BID)(x::Real)` is the same as `convert($($BID), x)`.

    `$($BID)(x::AbstractString)` is the same as `parse($($BID), x)`. This is provided for
    convenience since decimal literals are converted to `Float64` when parsed and may not
    produce what you expect.

    `$($BID)(sign, significand, exponent)` returns `sign * significand * 10^exponent`.
    If `sign` isn't passed, use the sign of `significand`.

    # Examples
    ```julia-repl
    julia> $($BID)(1)
    1.0

    julia> $($BID)(1.5)
    1.5

    julia> $($BID)("0.99999999999999999999999999999999999")
    1.0

    julia> $($BID)("0.99999999999999999999999999999999999", RoundDown)
    0.$("9"^(9*$w÷32-2))

    julia> $($BID)(-1, 123456, -4)
    -12.3456
    ```
    """ $BID

    @eval $BID(x::AbstractIrrational, r::RoundingMode) = $BID(string(BigFloat(x, precision=256)), r)

    @eval $BID(x::Union{Int8,UInt8,Int16,UInt16}) = $BID(Int32(x))
    @eval $BID(x::Float16) = $BID(Float32(x))

    @eval $BID(x::Rational{T}) where {T} = $BID(x.num) / $BID(x.den)

    # fallback conversions via strings
    @eval $BID(x::AbstractFloat) = $BID(string(x))
    @eval $BID(x::Integer) = $BID(string(x))

    @eval _precision(::Type{$BID}) = $(w == 32 ? 7 : w == 64 ? 16 : 34)
end

# match semantics in Julia 1.8
function Base.precision(::Union{T,Type{T}}; base::Integer=2) where {T<:DecimalFloatingPoint}
    base > 1 || throw(DomainError(base, "`base` cannot be less than 2."))
    p = _precision(T)
    return base == 10 ? Int(p) : floor(Int, p / log10(base))
end

# quickly check whether s begins with "±nan"
function isnanstr(s::AbstractString)
    st = iterate(s)
    c, i = '\0', 0
    while st !== nothing
        c, i = st
        isspace(c) || break
        st = iterate(s, i)
    end
    st === nothing && return false
    if c == '+' || c == '-'
        st = iterate(s, i)
        st === nothing && return false
        c, i = st
    end
    (c == 'n' || c == 'N') || return false
    st = iterate(s, i)
    if st !== nothing
        c, i = st
        if c == 'a' || c == 'A'
            st = iterate(s, i)
            if st !== nothing
                c, i = st
                if c == 'n' || c == 'N'
                    st = iterate(s, i)
                    if st === nothing
                        return true
                    end
                end
            end
        end
    end
    return false
end

"""
    sigexp(x::DecFP.DecimalFloatingPoint)

Return `(sign, significand, exponent)` such that `x` is equal to `sign * significand * 10^exponent`.
Throws `DomainError` for infinite or `NaN` arguments.

# Examples
```jldoctest
julia> sigexp(Dec64(1.25))
(1, 125, -2)
```
"""
sigexp(x::DecimalFloatingPoint)

"""
    exponent10(x::DecFP.DecimalFloatingPoint)

Get the exponent of the base 10 representation of a normalized floating-point number.

# Examples
```jldoctest
julia> exponent10(Dec64(123))
2
```
"""
exponent10(x::DecimalFloatingPoint)

"""
    ldexp10(x::DecFP.DecimalFloatingPoint, n::Integer)

Compute ``x * 10^n``.

# Examples
```jldoctest
julia> ldexp10(Dec64(15), 2)
1500.0
```
"""
ldexp10(x::DecimalFloatingPoint, n::Integer)

function Base.parse(::Type{T}, s::AbstractString) where {T<:DecimalFloatingPoint}
    x = _parse(T, s)
    if isnan(x) && !isnanstr(s)
        throw(ArgumentError("invalid number format $s"))
    end
    return x
end

function Base.tryparse(::Type{T}, s::AbstractString) where {T<:DecimalFloatingPoint}
    x = _parse(T, s)
    if isnan(x) && !isnanstr(s)
        return nothing
    end
    return x
end

Base.tryparse_internal(::Type{T}, s::AbstractString, startpos::Int, endpos::Int) where {T<:DecimalFloatingPoint} =
    tryparse(T, s[startpos:endpos]) # use a substring once Intel supports it?

_parse(::Type{T}, s::AbstractString) where {T<:DecimalFloatingPoint} =
    _parse(T, String(s)) # Intel only supports NUL-terminated strings

for w in (32,64,128)
    BID = Symbol(string("Dec",w))
    Ti = eval(Symbol(string("UInt",w)))
    Tsi = eval(Symbol(string("Int",w)))
    T = eval(BID)

    @eval _parse(::Type{$BID}, s::String) = ccall(($(bidsym(w,"from_string")), libbid), $BID, (Cstring,Cuint,Ref{Cuint}), s, _roundingmode(), zero(Cuint))

    @eval begin
        $BID(x::AbstractString) = parse($BID, x)

        function $BID(x::AbstractString, mode::RoundingMode)
            setrounding($BID, mode) do
                parse($BID, x)
            end
        end

        function tostring(x::$BID)
            # fills global _stringbuffer
            ccall(($(bidsym(w,"to_string")), libbid), Cvoid, (Ptr{UInt8},$BID,Ref{Cuint}), _stringbuffer(), x, zero(Cuint))
        end

        function Base.show(io::IO, x::$BID)
            isnan(x) && (print(io, "NaN"); return)
            isinf(x) && (print(io, signbit(x) ? "-Inf" : "Inf"); return)
            x == 0 && (print(io, signbit(x) ? "-0.0" : "0.0"); return)
            tostring(x)
            buffer = _stringbuffer()
            if buffer[1] == UInt8('-')
                print(io, '-')
            end
            normalized_exponent = exponent10(x)
            lastdigitindex = findfirst(isequal(UInt8('E')), buffer) - 1
            lastnonzeroindex = findlast(!isequal(UInt8('0')), view(buffer, 1:lastdigitindex))
            if -5 < normalized_exponent < 6
                # %f
                if normalized_exponent >= 0
                    if normalized_exponent >= lastnonzeroindex - 2
                        GC.@preserve buffer unsafe_write(io, pointer(buffer, 2), lastnonzeroindex - 1)
                        printzeros(io, normalized_exponent - lastnonzeroindex + 2)
                        print(io, ".0")
                    else
                        GC.@preserve buffer unsafe_write(io, pointer(buffer, 2), normalized_exponent + 1)
                        print(io, '.')
                        GC.@preserve buffer unsafe_write(io, pointer(buffer, normalized_exponent + 3), lastnonzeroindex - normalized_exponent - 2)
                    end
                else
                    print(io, "0.")
                    printzeros(io, -normalized_exponent - 1)
                    GC.@preserve buffer unsafe_write(io, pointer(buffer, 2), lastnonzeroindex - 1)
                end
            else
                # %e
                print(io, Char(buffer[2]), '.')
                if lastnonzeroindex == 2
                    print(io, '0')
                else
                    GC.@preserve buffer unsafe_write(io, pointer(buffer, 3), lastnonzeroindex - 2)
                end
                print(io, 'e')
                if normalized_exponent < 0
                    print(io, '-')
                    normalized_exponent = -normalized_exponent
                end
                b_lb = div(normalized_exponent, 10)
                b = 1
                while b <= b_lb
                    b *= 10
                end
                r = normalized_exponent
                while b > 0
                    q, r = divrem(r, b)
                    print(io, '0' + q)
                    b = div(b, 10)
                end
            end
            return
        end

        function sigexp(x::$BID)
            isnan(x) && throw(DomainError(x, "sigexp(x) is not defined for NaN."))
            isinf(x) && throw(DomainError(x, "sigexp(x) is only defined for finite x."))
            p = 9 * $w ÷ 32 - 2
            emax = 3 * 2^($w ÷ 16 + 3)
            bias = emax + p - 2
            ebits = $w ÷ 16 + 6
            sbits = $w - ebits - 1
            n = reinterpret($Ti, x)
            if n & $Ti(0x3) << ($w - 3) == $Ti(0x3) << ($w - 3)
                s = ((n & typemax($Ti) >> (ebits + 3)) | one($Ti) << sbits) % $Tsi
                e = ((n & typemax($Ti) << (sbits + 1) >> 3) >> (sbits - 2)) % Int
            else
                s = (n & typemax($Ti) >> (ebits + 1)) % $Tsi
                e = ((n & typemax($Ti) << (sbits + 1) >> 1) >> sbits) % Int
            end
            e -= bias
            return signbit(x) ? -1 : 1, s, e
        end

        Base.fma(x::$BID, y::$BID, z::$BID) = ccall(($(bidsym(w,"fma")), libbid), $BID, ($BID,$BID,$BID,Cuint,Ref{Cuint}), x, y, z, _roundingmode(), zero(Cuint))
        Base.muladd(x::$BID, y::$BID, z::$BID) = fma(x,y,z) # faster than x+y*z

        Base.one(::Union{Type{$BID},$BID}) = $(parse(T, "1"))
        Base.zero(::Union{Type{$BID},$BID}) = $(parse(T, "0"))

        Base.signbit(x::$BID) = $(zero(Ti)) != $(Ti(1) << (Ti(w - 1))) & reinterpret($Ti, x)
        Base.sign(x::$BID) = ifelse(isnan(x) || iszero(x), x, ifelse(signbit(x), $(parse(T, "-1")), $(parse(T, "1"))))

        Base.nextfloat(x::$BID) = ccall(($(bidsym(w,"nexttoward")), libbid), $BID, ($BID,Dec128,Ref{Cuint}), x, pinf128, zero(Cuint))
        Base.prevfloat(x::$BID) = ccall(($(bidsym(w,"nexttoward")), libbid), $BID, ($BID,Dec128,Ref{Cuint}), x, minf128, zero(Cuint))

        function Base.eps(x::$BID)
            if isfinite(x)
                # diff = nextfloat(x) - x, but check flags (and inline so that Ref is not heap-allocated)
                flags = Ref(zero(Cuint))
                next = ccall(($(bidsym(w,"nexttoward")), libbid), $BID, ($BID,Dec128,Ref{Cuint}), x, pinf128, flags)
                diff = ccall(($(bidsym(w,"sub")), libbid), $BID, ($BID,$BID,Cuint,Ref{Cuint}), next, x, _roundingmode(), flags)
                return @xchk(diff, flags[], OverflowError, "$($BID) value overflow", mask=OVERFLOW)
            else
                return $(parse(T, "NaN"))
            end
        end

        # the meaning of the exponent is different than for binary FP: it is 10^n, not 2^n:
        exponent10(x::$BID) = ccall(($(bidsym(w,"ilogb")), libbid), Cint, ($BID,Ref{Cuint}), x, zero(Cuint))
        ldexp10(x::$BID, n::Integer) = ccall(($(bidsym(w,"ldexp")), libbid), $BID, ($BID,Cint,Cuint,Ref{Cuint}), x, n, _roundingmode(), zero(Cuint))
    end

    for (f,c) in ((:isnan,"isNaN"), (:isinf,"isInf"), (:isfinite,"isFinite"), (:issubnormal,"isSubnormal"))
        @eval Base.$f(x::$BID) = ccall(($(bidsym(w,c)), libbid), Cint, ($BID,), x) != 0
    end

    for (f,c) in ((:+,"add"), (:-,"sub"), (:*,"mul"), (:/, "div"), (:hypot,"hypot"), (:atan,"atan2"), (:^,"pow"))
        @eval Base.$f(x::$BID, y::$BID) = ccall(($(bidsym(w,c)), libbid), $BID, ($BID,$BID,Cuint,Ref{Cuint}), x, y, _roundingmode(), zero(Cuint))
    end
    @eval Base.copysign(x::$BID, y::$BID) = ccall(($(bidsym(w,"copySign")), libbid), $BID, ($BID,$BID,Ref{Cuint}), x, y, zero(Cuint))

    for f in (:exp,:log,:sin,:cos,:tan,:asin,:acos,:atan,:sinh,:cosh,:tanh,:asinh,:acosh,:atanh,:log1p,:expm1,:log10,:log2,:exp2,:exp10,:sqrt,:cbrt)
        @eval Base.$f(x::$BID) = (flags = Ref(zero(Cuint)); @xchk(ccall(($(bidsym(w,f)), libbid), $BID, ($BID,Cuint,Ref{Cuint}), x, _roundingmode(), flags), flags[], DomainError, x, mask=INVALID))
    end
    @eval Base.abs(x::$BID) = (flags = Ref(zero(Cuint)); @xchk(ccall(($(bidsym(w,"abs")), libbid), $BID, ($BID,Ref{Cuint}), x, flags), flags[], DomainError, x, mask=INVALID))

    @eval Base.rem(x::$BID, y::$BID) = ccall(($(bidsym(w,"fmod")), libbid), $BID, ($BID,$BID,Ref{Cuint}), x, y, zero(Cuint))
    @eval Base.rem(x::$BID, y::$BID, ::RoundingMode{:Nearest}) = ccall(($(bidsym(w,"rem")), libbid), $BID, ($BID,$BID,Ref{Cuint}), x, y, zero(Cuint))

    @eval begin
        function Base.modf(x::$BID)
            ipart = Ref{$BID}()
            fpart = ccall(($(bidsym(w,"modf")), libbid), $BID, ($BID,Ref{$BID},Ref{Cuint}), x, ipart, zero(Cuint))
            fpart, ipart[]
        end
    end

    for (f,c) in ((:trunc,"round_integral_zero"), (:floor,"round_integral_negative"), (:ceil,"round_integral_positive"))
        @eval Base.$f(x::$BID) = (flags = Ref(zero(Cuint)); @xchk(ccall(($(bidsym(w,c)), libbid), $BID, ($BID,Ref{Cuint}), x, flags), flags[], DomainError, x, mask=INVALID))
    end
    @eval Base.:-(x::$BID) = ccall(($(bidsym(w,"negate")), libbid), $BID, ($BID,), x)
    @eval Base.round(x::$BID) = (flags = Ref(zero(Cuint)); @xchk(ccall(($(bidsym(w,"nearbyint")), libbid), $BID, ($BID,Cuint,Ref{Cuint}), x, _roundingmode(), flags), flags[], DomainError, x, mask=INVALID))

    @eval function SpecialFunctions.logabsgamma(x::$BID)
        isequal(modf(x)[1], -zero(x)) && return typemax(x), 1
        signgam = signbit(x) && mod(x, 2) > 1 ? -1 : 1
        flags = Ref(zero(Cuint))
        y = @xchk(ccall(($(bidsym(w,:lgamma)), libbid), $BID, ($BID,Cuint,Ref{Cuint}), x, _roundingmode(), flags), flags[], DomainError, x, mask=INVALID)
        return y, signgam
    end

    @eval SpecialFunctions.gamma(x::$BID) = (flags = Ref(zero(Cuint)); @xchk(ccall(($(bidsym(w,:tgamma)), libbid), $BID, ($BID,Cuint,Ref{Cuint}), x, _roundingmode(), flags), flags[], DomainError, x, mask=INVALID))

    @eval Random.rand(r::Random.AbstractRNG, ::Random.SamplerTrivial{Random.CloseOpen01{$BID}}) =  $BID(1, rand(r, zero($Ti):$Ti(maxintfloat($BID)) - one($Ti)), -(9 * $w ÷ 32 - 2))

    for (r,c) in ((RoundingMode{:Nearest},"round_integral_nearest_even"), (RoundingMode{:NearestTiesAway},"round_integral_nearest_away"), (RoundingMode{:ToZero},"round_integral_zero"), (RoundingMode{:Up},"round_integral_positive"), (RoundingMode{:Down},"round_integral_negative"))
        @eval Base.round(x::$BID, ::$r) = (flags = Ref(zero(Cuint)); @xchk(ccall(($(bidsym(w,c)), libbid), $BID, ($BID,Ref{Cuint}), x, flags), flags[], DomainError, x, mask=INVALID))
    end

    for (f,c) in ((:(==),"quiet_equal"), (:>,"quiet_greater"), (:<,"quiet_less"), (:(>=), "quiet_greater_equal"), (:(<=), "quiet_less_equal"))
        @eval Base.$f(x::$BID, y::$BID) = ccall(($(bidsym(w,c)), libbid), Cint, ($BID,$BID,Ref{Cuint}), x, y, zero(Cuint)) != 0
    end

    for Tf in (Float32,Float64)
        bT = string("binary",sizeof(Tf)*8)
        @eval begin
            Base.$(Symbol("$Tf"))(x::$BID) = ccall(($(bidsym(w,"to_",bT)), libbid), $Tf, ($BID,Cuint,Ref{Cuint}), x, _roundingmode(), zero(Cuint))
            $BID(x::$Tf) = ccall(($(string("__",bT,"_to_","bid",w)), libbid), $BID, ($Tf,Cuint,Ref{Cuint}), x, _roundingmode(), zero(Cuint))
        end
    end

    for c in (:π, :e, :ℯ, :γ, :catalan, :φ)
        @eval begin
            $BID(::Irrational{$(QuoteNode(c))}) = $(parse(T, setprecision(256) do
                                                                                      string(BigFloat(getfield(MathConstants, c)))
                                                                                  end))
        end
    end

    @eval promote_rule(::Type{$BID}, ::Type{Irrational}) = $BID

    for w′ in (32,64,128)
        BID′ = Symbol(string("Dec",w′))
        if w > w′
            @eval promote_rule(::Type{$BID}, ::Type{$BID′}) = $BID
        end
        if w > w′
            @eval $BID(x::$BID′) = ccall(($(string("__bid",w′,"_to_","bid",w)), libbid), $BID, ($BID′,Ref{Cuint}), x, zero(Cuint))
        elseif w < w′
            @eval $BID(x::$BID′) = ccall(($(string("__bid",w′,"_to_","bid",w)), libbid), $BID, ($BID′,Cuint,Ref{Cuint}), x, _roundingmode(), zero(Cuint))
        end

        # promote binary*decimal -> decimal, for consistency with other operations above
        # (there doesn't seem to be any clear standard for this)
        if w′ <= 64
            FP′ = Symbol(string("Float",w′))
            @eval promote_rule(::Type{$BID}, ::Type{$FP′}) = $(Symbol(string("Dec",max(w,w′))))
            for (i′, i′str) in (("Int$w′", "int$w′"), ("UInt$w′", "uint$w′"))
                Ti′ = eval(Symbol(i′))
                if w > w′
                    @eval $BID(x::$Ti′) = ccall(($(bidsym(w,"from_",i′str)), libbid), $BID, ($Ti′,Ref{Cuint}), x, zero(Cuint))
                else
                    @eval $BID(x::$Ti′) = ccall(($(bidsym(w,"from_",i′str)), libbid), $BID, ($Ti′,Cuint,Ref{Cuint}), x, _roundingmode(), zero(Cuint))
                end
            end
        end
    end

    for w′ in (8,16,32,64)
        for (i′, i′str) in (("Int$w′", "int$w′"), ("UInt$w′", "uint$w′"))
            Ti′ = eval(Symbol(i′))
            @eval begin
                Base.trunc(::Type{$Ti′}, x::$BID) = (flags = Ref(zero(Cuint)); @xchk(ccall(($(bidsym(w,"to_",i′str,"_xint")), libbid), $Ti′, ($BID,Ref{Cuint}), x, flags), flags[], InexactError, :trunc, $BID, x, mask=INVALID | OVERFLOW))
                Base.floor(::Type{$Ti′}, x::$BID) = (flags = Ref(zero(Cuint)); @xchk(ccall(($(bidsym(w,"to_",i′str,"_xfloor")), libbid), $Ti′, ($BID,Ref{Cuint}), x, flags), flags[], InexactError, :floor, $BID, x, mask=INVALID | OVERFLOW))
                Base.ceil(::Type{$Ti′}, x::$BID) = (flags = Ref(zero(Cuint)); @xchk(ccall(($(bidsym(w,"to_",i′str,"_xceil")), libbid), $Ti′, ($BID,Ref{Cuint}), x, flags), flags[], InexactError, :ceil, $BID, x, mask=INVALID | OVERFLOW))
                Base.round(::Type{$Ti′}, x::$BID, ::RoundingMode{:NearestTiesAway}) = (flags = Ref(zero(Cuint)); @xchk(ccall(($(bidsym(w,"to_",i′str,"_xrninta")), libbid), $Ti′, ($BID,Ref{Cuint}), x, flags), flags[], InexactError, :round, $BID, x, mask=INVALID | OVERFLOW))
                Base.$(Symbol("$Ti′"))(x::$BID) = (x::$BID) = (flags = Ref(zero(Cuint)); @xchk(ccall(($(bidsym(w,"to_",i′str,"_xfloor")), libbid), $Ti′, ($BID,Ref{Cuint}), x, flags), flags[], InexactError, :convert, $BID, x))
            end

            if w′ < w # integer conversion is exact
                @eval Base.:(==)(dec::$BID, int::$Ti′) = dec == $BID(int)
            end
        end
    end

    @eval Base.write(io::IO, x::$BID) = write(io, reinterpret($Ti, x))
    @eval Base.read(io::IO, x::Type{$BID}) = reinterpret($BID, read(io, $Ti))
    @eval Base.Float16(x::$BID) = Float16(Float32(x))
end # widths w

Base.round(x::DecimalFloatingPoint, ::RoundingMode{:FromZero}) = signbit(x) ? floor(x) : ceil(x)

for (f) in (:trunc, :floor, :ceil)
    @eval Base.$f(::Type{Signed}, x::DecimalFloatingPoint) = $f(Int, x)
    @eval Base.$f(::Type{Unsigned}, x::DecimalFloatingPoint) = $f(UInt, x)
    @eval Base.$f(::Type{Integer}, x::DecimalFloatingPoint) = $f(Int, x)

    @eval function Base.$f(::Type{I}, x::DecimalFloatingPoint) where {I<:Integer}
        x′ = $f(x)
        typemin(I) <= x′ <= typemax(I) || throw(InexactError(Symbol($f), I, x))
        _, s, e = sigexp(x′)
        return I(flipsign(s * I(10)^e, x))
    end
end

Base.Signed(x::DecimalFloatingPoint) = Int(x)
Base.Unsigned(x::DecimalFloatingPoint) = UInt(x)
Base.Integer(x::DecimalFloatingPoint) = Int(x)

function (::Type{I})(x::DecimalFloatingPoint) where {I<:Integer}
    x != trunc(x) && throw(InexactError(:convert, I, x))
    typemin(I) <= x <= typemax(I) || throw(InexactError(:convert, I, x))
    _, s, e = sigexp(x)
    return I(flipsign(s * I(10)^e, x))
end

Base.round(::Type{T}, x::DecimalFloatingPoint) where {T<:Integer} = convert(T, round(x))
Base.round(::Type{T}, x::DecimalFloatingPoint, ::RoundingMode{:Nearest}) where {T<:Integer} = convert(T, round(x, RoundNearest))
Base.round(::Type{T}, x::DecimalFloatingPoint, ::RoundingMode{:NearestTiesAway}) where {T<:Integer} = convert(T, round(x, RoundNearestTiesAway))
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

# see issue #92 — the fallback power_by_squaring fails for p < 0, and this is more accurate:
Base.:^(x::DecimalFloatingPoint, p::Integer) = x^oftype(x, p)
@inline Base.literal_pow(::typeof(^), x::DecimalFloatingPoint, ::Val{0}) = one(x)
@inline Base.literal_pow(::typeof(^), x::DecimalFloatingPoint, ::Val{1}) = x
@inline Base.literal_pow(::typeof(^), x::DecimalFloatingPoint, ::Val{2}) = x*x
@inline Base.literal_pow(::typeof(^), x::DecimalFloatingPoint, ::Val{3}) = x*x*x

# decompose is needed for hash() and comparison with Rational
function Base.decompose(x::DecimalFloatingPoint)::Tuple{BigInt, Int, BigInt}
    isnan(x) && return 0, 0, 0
    isinf(x) && return ifelse(signbit(x), -1, 1), 0, 0
    iszero(x) && return 0, 0, ifelse(signbit(x), -1, 1)
    sign, s, e = sigexp(x)
    if e >= 0
        if e <= 27
            return s * BigInt(Int64(5)^e), e, sign
        else
            return s * BigInt(5)^e, e, sign
        end
    else
        e2 = -e
        q, r = divrem(s, oftype(s, 5))
        while e2 > 0 && r == 0
            s = q
            q, r = divrem(s, oftype(s, 5))
            e2 -= 1
        end
        if e2 <= 27
            return copysign(s, x), e, Int64(5)^e2
        else
            return copysign(s, x), e, BigInt(5)^e2
        end
    end
end

function Base.:(==)(dec::DecimalFloatingPoint, rat::Rational)
    if isfinite(dec)
        rat.den == 1 && return dec == rat.num
        t = rat.den
        t >>= trailing_zeros(t)
        q, r = divrem(t, 5)
        while r == 0
            t = q
            q, r = divrem(t, 5)
        end
        return t == 1 && dec*rat.den == rat.num
    else
        return dec == rat.num/rat.den
    end
end

Base.:(==)(dec::T, num::Union{BigFloat,Float16,Float32,Float64,Integer}) where {T<:DecimalFloatingPoint} = dec == T(num, RoundUp) == T(num, RoundDown)
Base.:>(dec::T, num::Union{BigFloat,Float16,Float32,Float64,Integer,AbstractIrrational}) where {T<:DecimalFloatingPoint} = dec > T(num, RoundDown)
Base.:<(dec::T, num::Union{BigFloat,Float16,Float32,Float64,Integer,AbstractIrrational}) where {T<:DecimalFloatingPoint} = dec < T(num, RoundUp)
Base.:(>=)(dec::T, num::Union{BigFloat,Float16,Float32,Float64,Integer}) where {T<:DecimalFloatingPoint} = dec >= T(num, RoundUp)
Base.:(<=)(dec::T, num::Union{BigFloat,Float16,Float32,Float64,Integer}) where {T<:DecimalFloatingPoint} = dec <= T(num, RoundDown)

# canonicalize comparison order:
Base.:(==)(num::Union{BigFloat,Float16,Float32,Float64,Integer}, dec::T) where {T<:DecimalFloatingPoint} = dec == num
Base.:>(num::Union{BigFloat,Float16,Float32,Float64,Integer,AbstractIrrational}, dec::T) where {T<:DecimalFloatingPoint} = dec < num
Base.:<(num::Union{BigFloat,Float16,Float32,Float64,Integer,AbstractIrrational}, dec::T) where {T<:DecimalFloatingPoint} = dec > num
Base.:(>=)(num::Union{BigFloat,Float16,Float32,Float64,Integer}, dec::T) where {T<:DecimalFloatingPoint} = dec <= num
Base.:(<=)(num::Union{BigFloat,Float16,Float32,Float64,Integer}, dec::T) where {T<:DecimalFloatingPoint} = dec >= num

# used for next/prevfloat:
const pinf128 = parse(Dec128, "+Inf")
const minf128 = parse(Dec128, "-Inf")

for T in (Dec32, Dec64, Dec128)
    @eval begin
        Base.eps(::Type{$T}) = $(nextfloat(one(T)) - one(T))
        Base.typemax(::Type{$T}) = $(parse(T, "+inf"))
        Base.typemin(::Type{$T}) = $(parse(T, "-inf"))
        Base.floatmax(::Type{$T}) = $(prevfloat(parse(T, "+inf")))
    end
end

Base.floatmin(::Type{Dec32}) = reinterpret(Dec32, 0x03000001) # Dec32("1.0e-95")
Base.floatmin(::Type{Dec64}) = reinterpret(Dec64, 0x01e0000000000001) # Dec64("1.0e-383")
Base.floatmin(::Type{Dec128}) = reinterpret(Dec128, 0x00420000000000000000000000000001) # Dec128("1.0e-6143")

Base.maxintfloat(::Type{Dec32}) = reinterpret(Dec32, 0x36000001) # Dec32("1e7")
Base.maxintfloat(::Type{Dec64}) = reinterpret(Dec64, 0x33c0000000000001) # Dec64("1e16")
Base.maxintfloat(::Type{Dec128}) = reinterpret(Dec128, 0x30840000000000000000000000000001) # Dec128("1e34")
_int_maxintfloat(::Type{Dec32}) = 0x00989680
_int_maxintfloat(::Type{Dec64}) = 0x002386f26fc10000
_int_maxintfloat(::Type{Dec128}) = 0x0001ed09bead87c0378d8e6400000000

promote_rule(::Type{F}, ::Type{Float16}) where {F<:DecimalFloatingPoint} = F
promote_rule(::Type{F}, ::Type{T}) where {F<:DecimalFloatingPoint,T<:Union{Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64,Int128,UInt128}} = F

Base.BigFloat(x::DecimalFloatingPoint) = BigFloat(string(x))
Base.BigFloat(x::DecimalFloatingPoint, r::RoundingMode; precision::Integer=Base.MPFR.DEFAULT_PRECISION[]) = BigFloat(string(x), r; precision=precision)

# so that mathconsts get promoted to Dec32, not Dec64, like Float32
promote_rule(::Type{Irrational{s}}, ::Type{F}) where {s,F<:DecimalFloatingPoint} = F
promote_rule(::Type{Irrational{s}}, T::Type{Complex{F}}) where {s,F<:DecimalFloatingPoint} = T

Base.widen(::Type{Dec32}) = Dec64
Base.widen(::Type{Dec64}) = Dec128

include("printf.jl")

macro d_str(s, flags...) parse(Dec64, s) end
macro d32_str(s, flags...) parse(Dec32, s) end
macro d64_str(s, flags...) parse(Dec64, s) end
macro d128_str(s, flags...) parse(Dec128, s) end

# for zero-padding in printing routines above
function printzeros(io::IO, n::Int)
    for i = 1:n
        print(io, '0')
    end
end

end # module
