
module DecFP

using DecFP_jll

import Printf, SpecialFunctions

export Dec32, Dec64, Dec128, @d_str, @d32_str, @d64_str, @d128_str, exponent10, ldexp10

const _buffer = Vector{Vector{UInt8}}(undef, Threads.nthreads())

import Base.promote_rule

# flags isn't defined until __init__ runs
const _flags = [0x00000000]
const flags = Vector{Vector{Cuint}}(undef, Threads.nthreads())

# clear exception flags and return x
function nox(x)
    flags[Threads.threadid()][1] = 0
    return x
end

# Check exception flags in mask & throw, otherwise returning x;
# always clearing exceptions.  This is a macros so that
# the error message is only evaluated if an exception occurs.
macro xchk(x, exc, args...)
    mask=0x3f
    if !isempty(args) && Meta.isexpr(args[end], :(=)) && args[end].args[1] == :mask # mask=... keyword at end
        mask = esc(args[end].args[2])
        args = args[1:end-1]
    end
    quote
        ret = $(esc(x))
        if $exc === nothing
            flags[Threads.threadid()][1] = 0
        else
            f = flags[Threads.threadid()][1]
            flags[Threads.threadid()][1] = 0
            f & $mask != 0 && throw($exc($(map(esc,args)...)))
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

const roundingmode = Vector{DecFPRoundingMode}(undef, Threads.nthreads())

# global vectors must be initialized at runtime (via __init__)
function __init__()
    for i = 1:Threads.nthreads()
        global _buffer[i] = fill(0x00, 1024)
        global flags[i] = [0x00000000]
        global roundingmode[i] = DecFPRoundNearest
    end
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

Base.Rounding.rounding_raw(::Type{T}) where {T<:DecimalFloatingPoint} =
    roundingmode[Threads.threadid()]
Base.Rounding.setrounding_raw(::Type{T}, r::DecFPRoundingMode) where {T<:DecimalFloatingPoint} =
    roundingmode[Threads.threadid()] = r

Base.Rounding.rounding(::Type{T}) where {T<:DecimalFloatingPoint} =
    convert(RoundingMode, Base.Rounding.rounding_raw(T))
Base.Rounding.setrounding(::Type{T}, r::RoundingMode) where {T<:DecimalFloatingPoint} =
    Base.Rounding.setrounding_raw(T, convert(DecFPRoundingMode, r))

for w in (32,64,128)
    BID = Symbol(string("Dec",w))
    Ti = Symbol(string("UInt",w))
    @eval struct $BID <: DecimalFloatingPoint
        x::$Ti
        $BID(x::Number) = convert($BID, x)
        Base.reinterpret(::Type{$BID}, x::$Ti) = new(x)
    end
    # fix method ambiguities:
    @eval $BID(x::Rational{T}) where {T} = convert($BID, x)
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
ldexp10(x::DecFP.DecimalFloatingPoint, n::Integer)

# for compatibility with julia#29885
if isdefined(Base.Grisu, :getbuf)
    getdigitsbuf() = Base.Grisu.getbuf()
elseif isdefined(Base.Grisu, :DIGITSs)
    getdigitsbuf() = Base.Grisu.DIGITSs[Threads.threadid()]
else
    getdigitsbuf() = Base.Grisu.DIGITS
end

for w in (32,64,128)
    BID = Symbol(string("Dec",w))
    Ti = eval(Symbol(string("UInt",w)))
    T = eval(BID)

    # hack: we need an internal parsing function that doesn't check exceptions, since
    # flags isn't defined until __init__ runs.  Similarly for nextfloat/prevfloat
    @eval begin
        _parse(::Type{$BID}, s::AbstractString) =
            ccall(($(bidsym(w,"from_string")), libbid), $BID, (Ptr{UInt8},Cuint,Ref{Cuint}), s, DecFPRoundNearest, _flags)
        _nextfloat(x::$BID) = ccall(($(bidsym(w,"nexttoward")), libbid), $BID, ($BID,Dec128,Ref{Cuint}), x, pinf128, _flags)
        _prevfloat(x::$BID) = ccall(($(bidsym(w,"nexttoward")), libbid), $BID, ($BID,Dec128,Ref{Cuint}), x, minf128, _flags)
        _sub(x::$BID, y::$BID) = ccall(($(bidsym(w,"sub")), libbid), $BID, ($BID,$BID,Cuint,Ref{Cuint}), x, y, DecFPRoundNearest, _flags)
    end

    @eval begin
        function Base.parse(::Type{$BID}, s::AbstractString)
            x = _parse($BID, s)
            if isnan(x) && !isnanstr(s)
                throw(ArgumentError("invalid number format $s"))
            end
            return @xchk(x, nothing)
        end

        $BID(x::AbstractString) = parse($BID, x)

        function tostring(x::$BID)
            # fills global _buffer
            ccall(($(bidsym(w,"to_string")), libbid), Cvoid, (Ptr{UInt8},$BID,Ref{Cuint}), _buffer[Threads.threadid()], x, flags[Threads.threadid()])
        end

        function Base.show(io::IO, x::$BID)
            isnan(x) && (print(io, "NaN"); return)
            isinf(x) && (print(io, signbit(x) ? "-Inf" : "Inf"); return)
            x == 0 && (print(io, signbit(x) ? "-0.0" : "0.0"); return)
            tostring(x)
            if _buffer[Threads.threadid()][1] == UInt8('-')
                print(io, '-')
            end
            normalized_exponent = exponent10(x)
            lastdigitindex = findfirst(isequal(UInt8('E')), _buffer[Threads.threadid()]) - 1
            lastnonzeroindex = findlast(!isequal(UInt8('0')), view(_buffer[Threads.threadid()], 1:lastdigitindex))
            if -5 < normalized_exponent < 6
                # %f
                if normalized_exponent >= 0
                    if normalized_exponent >= lastnonzeroindex - 2
                        GC.@preserve _buffer unsafe_write(io, pointer(_buffer[Threads.threadid()], 2), lastnonzeroindex - 1)
                        printzeros(io, normalized_exponent - lastnonzeroindex + 2)
                        print(io, ".0")
                    else
                        GC.@preserve _buffer unsafe_write(io, pointer(_buffer[Threads.threadid()], 2), normalized_exponent + 1)
                        print(io, '.')
                        GC.@preserve _buffer unsafe_write(io, pointer(_buffer[Threads.threadid()], normalized_exponent + 3), lastnonzeroindex - normalized_exponent - 2)
                    end
                else
                    print(io, "0.")
                    printzeros(io, -normalized_exponent - 1)
                    GC.@preserve _buffer unsafe_write(io, pointer(_buffer[Threads.threadid()], 2), lastnonzeroindex - 1)
                end
            else
                # %e
                print(io, Char(_buffer[Threads.threadid()][2]), '.')
                if lastnonzeroindex == 2
                    print(io, '0')
                else
                    GC.@preserve _buffer unsafe_write(io, pointer(_buffer[Threads.threadid()], 3), lastnonzeroindex - 2)
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

        function Printf.Printf.fix_dec(x::$BID, n::Int, digits)
            if n > length(digits) - 1
                n = length(digits) - 1
            end
            rounded = round(ldexp10(x, n), RoundNearestTiesAway)
            if rounded == 0
                digits[1] = UInt8('0')
                return Int32(1), Int32(1), signbit(x)
            end
            tostring(rounded)
            trailing_zeros = 0
            i = 2
            while _buffer[Threads.threadid()][i] != UInt8('E')
                digits[i - 1] = _buffer[Threads.threadid()][i]
                if _buffer[Threads.threadid()][i] == UInt8('0')
                    trailing_zeros += 1
                else
                    trailing_zeros = 0
                end
                i += 1
            end
            ndigits = i - 2
            len = ndigits - trailing_zeros
            i += 1
            if _buffer[Threads.threadid()][i] == UInt8('+')
                expsign = +1
            elseif _buffer[Threads.threadid()][i] == UInt8('-')
                expsign = -1
            end
            exponent = 0
            i += 1
            while _buffer[Threads.threadid()][i] != 0x00
                exponent = exponent * 10 + _buffer[Threads.threadid()][i] - UInt8('0')
                i += 1
            end
            exponent *= expsign
            pt = ndigits + exponent - n
            neg = signbit(x)
            return Int32(len), Int32(pt), neg
        end

        function Printf.Printf.ini_dec(x::$BID, n::Int, digits)
            if n > length(digits) - 1
                n = length(digits) - 1
            end
            if x == 0
                for i = 1:n
                    digits[i] = UInt8('0')
                end
                return Int32(1), Int32(1), signbit(x)
            end
            normalized_exponent = exponent10(x)
            rounded = round(ldexp10(x, n - 1 - normalized_exponent), RoundNearestTiesAway)
            rounded_exponent = exponent10(rounded)
            tostring(rounded)
            i = 2
            while _buffer[Threads.threadid()][i] != UInt8('E')
                digits[i - 1] = _buffer[Threads.threadid()][i]
                i += 1
            end
            while i <= n + 1
                digits[i - 1] = UInt8('0')
                i += 1
            end
            pt = normalized_exponent + rounded_exponent - n + 2
            neg = signbit(x)
            return Int32(n), Int32(pt), neg
        end

        # compatibility with julia#30373
        Printf.Printf.fix_dec(x::$BID, n::Int) = Printf.Printf.fix_dec(x, n, getdigitsbuf())
        Printf.Printf.ini_dec(x::$BID, n::Int) = Printf.Printf.ini_dec(x, n, getdigitsbuf())

        Base.fma(x::$BID, y::$BID, z::$BID) = nox(ccall(($(bidsym(w,"fma")), libbid), $BID, ($BID,$BID,$BID,Cuint,Ref{Cuint}), x, y, z, roundingmode[Threads.threadid()], flags[Threads.threadid()]))
        Base.muladd(x::$BID, y::$BID, z::$BID) = fma(x,y,z) # faster than x+y*z

        Base.one(::Union{Type{$BID},$BID}) = $(_parse(T, "1"))
        Base.zero(::Union{Type{$BID},$BID}) = $(_parse(T, "0"))

        Base.signbit(x::$BID) = $(zero(Ti)) != $(Ti(1) << (Ti(w - 1))) & x.x
        Base.sign(x::$BID) = ifelse(signbit(x), $(_parse(T, "-1")), $(_parse(T, "1")))

        Base.nextfloat(x::$BID) = nox(_nextfloat(x))
        Base.prevfloat(x::$BID) = nox(_prevfloat(x))
        Base.eps(x::$BID) = ifelse(isfinite(x), @xchk(nextfloat(x) - x, OverflowError, "$($BID) value overflow", mask=OVERFLOW), $(_parse(T, "NaN")))

        # the meaning of the exponent is different than for binary FP: it is 10^n, not 2^n:
        exponent10(x::$BID) = nox(ccall(($(bidsym(w,"ilogb")), libbid), Cint, ($BID,Ref{Cuint}), x, flags[Threads.threadid()]))
        ldexp10(x::$BID, n::Integer) = nox(ccall(($(bidsym(w,"ldexp")), libbid), $BID, ($BID,Cint,Cuint,Ref{Cuint}), x, n, roundingmode[Threads.threadid()], flags[Threads.threadid()]))
    end

    for (f,c) in ((:isnan,"isNaN"), (:isinf,"isInf"), (:isfinite,"isFinite"), (:issubnormal,"isSubnormal"))
        @eval Base.$f(x::$BID) = ccall(($(bidsym(w,c)), libbid), Cint, ($BID,), x) != 0
    end

    for (f,c) in ((:+,"add"), (:-,"sub"), (:*,"mul"), (:/, "div"), (:hypot,"hypot"), (:atan,"atan2"), (:^,"pow"))
        @eval Base.$f(x::$BID, y::$BID) = nox(ccall(($(bidsym(w,c)), libbid), $BID, ($BID,$BID,Cuint,Ref{Cuint}), x, y, roundingmode[Threads.threadid()], flags[Threads.threadid()]))
    end
    @eval Base.copysign(x::$BID, y::$BID) = nox(ccall(($(bidsym(w,"copySign")), libbid), $BID, ($BID,$BID,Ref{Cuint}), x, y, flags[Threads.threadid()]))

    for f in (:exp,:log,:sin,:cos,:tan,:asin,:acos,:atan,:sinh,:cosh,:tanh,:asinh,:acosh,:atanh,:log1p,:expm1,:log10,:log2,:exp2,:exp10,:sqrt,:cbrt)
        @eval Base.$f(x::$BID) = @xchk(ccall(($(bidsym(w,f)), libbid), $BID, ($BID,Cuint,Ref{Cuint}), x, roundingmode[Threads.threadid()], flags[Threads.threadid()]), DomainError, x, mask=INVALID)
    end
    @eval Base.abs(x::$BID) = @xchk(ccall(($(bidsym(w,"abs")), libbid), $BID, ($BID,Ref{Cuint}), x, flags[Threads.threadid()]), DomainError, x, mask=INVALID)

    for (f,c) in ((:trunc,"round_integral_zero"), (:floor,"round_integral_negative"), (:ceil,"round_integral_positive"))
        @eval Base.$f(x::$BID) = @xchk(ccall(($(bidsym(w,c)), libbid), $BID, ($BID,Ref{Cuint}), x, flags[Threads.threadid()]), DomainError, x, mask=INVALID)
    end
    @eval Base.:-(x::$BID) = @xchk(ccall(($(bidsym(w,"negate")), libbid), $BID, ($BID,), x), DomainError, x, mask=INVALID)
    @eval Base.round(x::$BID) = @xchk(ccall(($(bidsym(w,"nearbyint")), libbid), $BID, ($BID,Cuint,Ref{Cuint}), x, roundingmode[Threads.threadid()], flags[Threads.threadid()]), DomainError, x, mask=INVALID)

    @eval SpecialFunctions.lgamma(x::$BID) = @xchk(ccall(($(bidsym(w,:lgamma)), libbid), $BID, ($BID,Cuint,Ref{Cuint}), x, roundingmode[Threads.threadid()], flags[Threads.threadid()]), DomainError, x, mask=INVALID)
    @eval SpecialFunctions.gamma(x::$BID) = @xchk(ccall(($(bidsym(w,:tgamma)), libbid), $BID, ($BID,Cuint,Ref{Cuint}), x, roundingmode[Threads.threadid()], flags[Threads.threadid()]), DomainError, x, mask=INVALID)

    for (r,c) in ((RoundingMode{:Nearest},"round_integral_nearest_even"), (RoundingMode{:NearestTiesAway},"round_integral_nearest_away"), (RoundingMode{:ToZero},"round_integral_zero"), (RoundingMode{:Up},"round_integral_positive"), (RoundingMode{:Down},"round_integral_negative"))
        @eval Base.round(x::$BID, ::$r) = @xchk(ccall(($(bidsym(w,c)), libbid), $BID, ($BID,Ref{Cuint}), x, flags[Threads.threadid()]), DomainError, x, mask=INVALID)
    end

    for (f,c) in ((:(==),"quiet_equal"), (:>,"quiet_greater"), (:<,"quiet_less"), (:(>=), "quiet_greater_equal"), (:(<=), "quiet_less_equal"))
        @eval Base.$f(x::$BID, y::$BID) = nox(ccall(($(bidsym(w,c)), libbid), Cint, ($BID,$BID,Ref{Cuint}), x, y, flags[Threads.threadid()]) != 0)
    end

    for Tf in (Float32,Float64)
        bT = string("binary",sizeof(Tf)*8)
        @eval begin
            Base.convert(::Type{$Tf}, x::$BID) = nox(ccall(($(bidsym(w,"to_",bT)), libbid), $Tf, ($BID,Cuint,Ref{Cuint}), x, roundingmode[Threads.threadid()], flags[Threads.threadid()]))
            Base.$(Symbol("$Tf"))(x::$BID) = convert($Tf, x)
            Base.convert(::Type{$BID}, x::$Tf) = nox(ccall(($(string("__",bT,"_to_","bid",w)), libbid), $BID, ($Tf,Cuint,Ref{Cuint}), x, roundingmode[Threads.threadid()], flags[Threads.threadid()]))
        end
    end

    for c in (:π, :e, :ℯ, :γ, :catalan, :φ)
        @eval begin
            Base.convert(::Type{$BID}, ::Irrational{$(QuoteNode(c))}) = $(_parse(T, setprecision(256) do
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
            @eval Base.convert(::Type{$BID}, x::$BID′) = @xchk(ccall(($(string("__bid",w′,"_to_","bid",w)), libbid), $BID, ($BID′,Ref{Cuint}), x, flags[Threads.threadid()]), nothing)
        elseif w < w′
            @eval Base.convert(::Type{$BID}, x::$BID′) = @xchk(ccall(($(string("__bid",w′,"_to_","bid",w)), libbid), $BID, ($BID′,Cuint,Ref{Cuint}), x, roundingmode[Threads.threadid()], flags[Threads.threadid()]), nothing)
        end

        # promote binary*decimal -> decimal, for consistency with other operations above
        # (there doesn't seem to be any clear standard for this)
        if w′ <= 64
            FP′ = Symbol(string("Float",w′))
            @eval promote_rule(::Type{$BID}, ::Type{$FP′}) = $(Symbol(string("Dec",max(w,w′))))
            for (i′, i′str) in (("Int$w′", "int$w′"), ("UInt$w′", "uint$w′"))
                Ti′ = eval(Symbol(i′))
                if w > w′
                    @eval Base.convert(::Type{$BID}, x::$Ti′) = nox(ccall(($(bidsym(w,"from_",i′str)), libbid), $BID, ($Ti′,Ref{Cuint}), x, flags[Threads.threadid()]))
                else
                    @eval Base.convert(::Type{$BID}, x::$Ti′) = nox(ccall(($(bidsym(w,"from_",i′str)), libbid), $BID, ($Ti′,Cuint,Ref{Cuint}), x, roundingmode[Threads.threadid()], flags[Threads.threadid()]))
                end
            end
        end
    end

    for w′ in (8,16,32,64)
        for (i′, i′str) in (("Int$w′", "int$w′"), ("UInt$w′", "uint$w′"))
            Ti′ = eval(Symbol(i′))
            @eval begin
                Base.trunc(::Type{$Ti′}, x::$BID) = @xchk(ccall(($(bidsym(w,"to_",i′str,"_xint")), libbid), $Ti′, ($BID,Ref{Cuint}), x, flags[Threads.threadid()]), InexactError, :trunc, $BID, x, mask=INVALID | OVERFLOW)
                Base.floor(::Type{$Ti′}, x::$BID) = @xchk(ccall(($(bidsym(w,"to_",i′str,"_xfloor")), libbid), $Ti′, ($BID,Ref{Cuint}), x, flags[Threads.threadid()]), InexactError, :floor, $BID, x, mask=INVALID | OVERFLOW)
                Base.ceil(::Type{$Ti′}, x::$BID) = @xchk(ccall(($(bidsym(w,"to_",i′str,"_xceil")), libbid), $Ti′, ($BID,Ref{Cuint}), x, flags[Threads.threadid()]), InexactError, :ceil, $BID, x, mask=INVALID | OVERFLOW)
                Base.round(::Type{$Ti′}, x::$BID, ::RoundingMode{:NearestTiesAway}) = @xchk(ccall(($(bidsym(w,"to_",i′str,"_xrninta")), libbid), $Ti′, ($BID,Ref{Cuint}), x, flags[Threads.threadid()]), InexactError, :round, $BID, x, mask=INVALID | OVERFLOW)
                Base.convert(::Type{$Ti′}, x::$BID) = @xchk(ccall(($(bidsym(w,"to_",i′str,"_xfloor")), libbid), $Ti′, ($BID,Ref{Cuint}), x, flags[Threads.threadid()]), InexactError, :convert, $BID, x)
                Base.$(Symbol("$Ti′"))(x::$BID) = convert($Ti′, x)
            end
        end
    end

    @eval Base.bswap(x::$BID) = reinterpret($BID, bswap(x.x))
    @eval Base.convert(::Type{Float16}, x::$BID) = convert(Float16, convert(Float32, x))
    @eval Base.Float16(x::$BID) = convert(Float16, x)
    @eval Base.reinterpret(::Type{$Ti}, x::$BID) = x.x
end # widths w

Base.round(x::DecimalFloatingPoint, ::RoundingMode{:FromZero}) = signbit(x) ? floor(x) : ceil(x)

Base.trunc(::Type{Integer}, x::DecimalFloatingPoint) = trunc(Int, x)
Base.floor(::Type{Integer}, x::DecimalFloatingPoint) = floor(Int, x)
Base.ceil(::Type{Integer}, x::DecimalFloatingPoint) = ceil(Int, x)
Base.round(::Type{Integer}, x::DecimalFloatingPoint) = round(Int, x)
Base.round(::Type{Integer}, x::DecimalFloatingPoint, ::RoundingMode{:NearestTiesAway}) = round(Int, x, RoundNearestTiesAway)
Base.convert(::Type{Integer}, x::DecimalFloatingPoint) = convert(Int, x)

Base.round(::Type{T}, x::DecimalFloatingPoint) where {T<:Integer} = convert(T, round(x))
Base.round(::Type{T}, x::DecimalFloatingPoint, ::RoundingMode{:Nearest}) where {T<:Integer} = convert(T, round(x, RoundNearest))
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

# used for next/prevfloat:
const pinf128 = _parse(Dec128, "+Inf")
const minf128 = _parse(Dec128, "-Inf")

for T in (Dec32, Dec64, Dec128)
    @eval begin
        Base.eps(::Type{$T}) = $(_sub(_nextfloat(one(T)), one(T)))
        Base.typemax(::Type{$T}) = $(_parse(T, "+inf"))
        Base.typemin(::Type{$T}) = $(_parse(T, "-inf"))
        Base.floatmax(::Type{$T}) = $(_prevfloat(_parse(T, "+inf")))
        Base.floatmin(::Type{$T}) = $(_nextfloat(zero(T)))
    end
end

Base.maxintfloat(::Type{Dec32}) = reinterpret(Dec32, 0x36000001) # Dec32("1e7")
Base.maxintfloat(::Type{Dec64}) = reinterpret(Dec64, 0x33c0000000000001) # Dec64("1e16")
Base.maxintfloat(::Type{Dec128}) = reinterpret(Dec128, 0x30840000000000000000000000000001) # Dec128("1e34")

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

# for zero-padding in printing routines above
function printzeros(io::IO, n::Int)
    for i = 1:n
        print(io, '0')
    end
end

end # module
