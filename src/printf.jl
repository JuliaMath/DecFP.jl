# max exponent: `ceil(Int, log10(floatmax(T)))`
max_integer_part_width(::Type{Dec32}) = 97
max_integer_part_width(::Type{Dec64}) = 385
max_integer_part_width(::Type{Dec128}) = 6145

Printf.plength(f::Printf.Spec{S}, x::T) where {S<:Printf.Floats, T<:DecimalFloatingPoint} = max(f.width, max_integer_part_width(T) + f.precision + 2)

Printf.tofloat(x::DecimalFloatingPoint) = x

@inline function Printf.fmt(buf, pos, x::DecimalFloatingPoint, spec::Printf.Spec{T}) where {T<:Printf.Floats}
    leftalign, plus, space, zero, hash, width, prec =
        spec.leftalign, spec.plus, spec.space, spec.zero, spec.hash, spec.width, spec.precision
    if T == Val{'e'} || T == Val{'E'}
        newpos = Base.Ryu.writeexp(buf, pos, x, prec, plus, space, hash, Printf.char(T), UInt8('.'))
    elseif T == Val{'f'} || T == Val{'F'}
        newpos = Base.Ryu.writefixed(buf, pos, x, prec, plus, space, hash, UInt8('.'))
    elseif T == Val{'g'} || T == Val{'G'}
        if isinf(x) || isnan(x)
            newpos = Base.Ryu.writeshortest(buf, pos, x, plus, space)
        else
            # C11-compliant general format
            prec = prec == 0 ? 1 : prec
            _, sig, dexp = sigexp(x)
            numdig = ndigits(sig)
            exp = dexp + numdig - 1
            if -4 ≤ exp < prec
                newpos = Base.Ryu.writefixed(buf, pos, x, prec - (exp + 1), plus, space, hash, UInt8('.'), !hash)
            else
                newpos = Base.Ryu.writeexp(buf, pos, x, prec - 1, plus, space, hash, T == Val{'g'} ? UInt8('e') : UInt8('E'), UInt8('.'), !hash)
            end
        end
    elseif T == Val{'a'} || T == Val{'A'}
        throw(ArgumentError("%a format specifier is not implemented for DecimalFloatingPoint"))
    end
    if newpos - pos < width
        # need to pad
        if leftalign
            # easy case, just pad spaces after number
            for _ = 1:(width - (newpos - pos))
                buf[newpos] = UInt8(' ')
                newpos += 1
            end
        else
            # right aligned
            n = width - (newpos - pos)
            if zero && isfinite(x)
                ex = (x < 0 || (plus | space)) + (T <: Union{Val{'a'}, Val{'A'}} ? 2 : 0)
                so = pos + ex
                len = (newpos - pos) - ex
                copyto!(buf, so + n, buf, so, len)
                for i = so:(so + n - 1)
                    buf[i] = UInt8('0')
                end
                newpos += n
            else
                copyto!(buf, pos + n, buf, pos, newpos - pos)
                for i = pos:(pos + n - 1)
                    buf[i] = UInt8(' ')
                end
                newpos += n
            end
        end
    end
    return newpos
end
        
function Base.Ryu.writefixed(buf, pos, x::T, precision=-1, plus=false, space=false, hash=false, decchar=UInt8('.'), trimtrailingzeros=false) where {T<:DecimalFloatingPoint}
    pos = Base.Ryu.append_sign(x, plus, space, buf, pos)

    # special cases
    if iszero(x)
        buf[pos] = UInt8('0')
        pos += 1
        if precision > 0 && !trimtrailingzeros
            buf[pos] = decchar
            pos += 1
            for _ = 1:precision
                buf[pos] = UInt8('0')
                pos += 1
            end
        elseif hash
            buf[pos] = decchar
            pos += 1
        end
        return pos
    elseif isnan(x)
        buf[pos] = UInt8('N')
        buf[pos + 1] = UInt8('a')
        buf[pos + 2] = UInt8('N')
        return pos + 3
    elseif !isfinite(x)
        buf[pos] = UInt8('I')
        buf[pos + 1] = UInt8('n')
        buf[pos + 2] = UInt8('f')
        return pos + 3
    end

    _, sig, exp = sigexp(x)
    numdig = ndigits(sig)
    digits = codeunits(string(sig))
    i = numdig
    while digits[i] == UInt8('0')
        sig ÷= 10
        exp += 1
        numdig -= 1
        i -= 1
    end
    if trimtrailingzeros && exp <= 0 && -exp < precision
        precision = -exp
    end
    reqdig = numdig + exp + precision
    if reqdig < 0
        buf[pos] = UInt8('0')
        pos += 1
        if precision > 0 && !trimtrailingzeros
            buf[pos] = decchar
            pos += 1
            for _ = 1:precision
                buf[pos] = UInt8('0')
                pos += 1
            end
        elseif hash
            buf[pos] = decchar
            pos += 1
        end
        return pos
    end
    if reqdig < numdig
        diff = numdig - reqdig
        if diff == 1
            denominator = oftype(sig, 10)
        elseif diff == 2
            denominator = oftype(sig, 100)
        elseif diff == 3
            denominator = oftype(sig, 1000)
        else
            denominator = oftype(sig, 10)^Int32(diff)
        end
        sig = div(sig, denominator, RoundNearest)
        exp += diff
    end
    digits = codeunits(string(sig))
    if exp <= -length(digits)
        buf[pos] = UInt8('0')
        buf[pos + 1] = decchar
        pos += 2
        for i in 1:-(exp + length(digits))
            buf[pos] = UInt8('0')
            pos += 1
        end
    end
    decindex = length(digits) + exp
    for i in 1:length(digits)
        buf[pos] = digits[i]
        pos += 1
        if i == decindex && (precision > 0 || hash)
            buf[pos] = decchar
            pos += 1
        end
    end
    for i in numdig+1:reqdig
        i > decindex && trimtrailingzeros && break
        buf[pos] = UInt8('0')
        pos += 1
        if i == decindex && (hash || i != reqdig)
            trimtrailingzeros && break
            buf[pos] = decchar
            pos += 1
        end
    end
    return pos
end

function Base.Ryu.writeexp(buf, pos, x::T, precision=-1, plus=false, space=false, hash=false, expchar=UInt8('e'), decchar=UInt8('.'), trimtrailingzeros=false) where {T<:DecimalFloatingPoint}
    pos = Base.Ryu.append_sign(x, plus, space, buf, pos)

    # special cases
    if iszero(x)
        @inbounds buf[pos] = UInt8('0')
        pos += 1
        if precision > 0 && !trimtrailingzeros
            @inbounds buf[pos] = decchar
            pos += 1
            for _ = 1:precision
                @inbounds buf[pos] = UInt8('0')
                pos += 1
            end
        elseif hash
            @inbounds buf[pos] = decchar
            pos += 1
        end
        @inbounds buf[pos] = expchar
        @inbounds buf[pos + 1] = UInt8('+')
        @inbounds buf[pos + 2] = UInt8('0')
        @inbounds buf[pos + 3] = UInt8('0')
        return pos + 4
    elseif isnan(x)
        @inbounds buf[pos] = UInt8('N')
        @inbounds buf[pos + 1] = UInt8('a')
        @inbounds buf[pos + 2] = UInt8('N')
        return pos + 3
    elseif !isfinite(x)
        @inbounds buf[pos] = UInt8('I')
        @inbounds buf[pos + 1] = UInt8('n')
        @inbounds buf[pos + 2] = UInt8('f')
        return pos + 3
    end

    _, sig, exp = sigexp(x)
    numdig = ndigits(sig)
    if numdig > precision + 1
        diff = numdig - precision - 1
        denominator = oftype(sig, 10)^Int32(diff)
        sig = div(sig, denominator, RoundNearest)
        exp += diff
    end
    numdig = ndigits(sig)
    if numdig > precision + 1
        # leading digit 9 rounded up to 10
        sig ÷= 10
        exp += 1
    end
    digits = codeunits(string(sig))
    numdig = length(digits)
    if trimtrailingzeros
        i = numdig
        while digits[i] == UInt8('0')
            exp += 1
            numdig -= 1
            i -= 1
        end
    end
    pexp = exp + numdig - 1
    buf[pos] = digits[1]
    pos += 1
    if numdig > 1 || (precision > 0 && !trimtrailingzeros) || hash
        buf[pos] = decchar
        pos += 1
    end
    for i in 2:numdig
        buf[pos] = digits[i]
        pos += 1
    end
    if !trimtrailingzeros
        for _ in numdig:precision
            buf[pos] = UInt8('0')
            pos += 1
        end
    end
    buf[pos] = expchar
    pos += 1
    if pexp < 0
        @inbounds buf[pos] = UInt8('-')
        pos += 1
        pexp = -pexp
    else
        @inbounds buf[pos] = UInt8('+')
        pos += 1
    end
    if pexp >= 1000
        c = (pexp % 100) % UInt8
        @inbounds d100 = Base._dec_d100[div(pexp, 100) + 1]
        @inbounds buf[pos] = d100 % UInt8
        @inbounds buf[pos + 1] = (d100 >> 0x8) % UInt8
        @inbounds d100 = Base._dec_d100[c + 1]
        @inbounds buf[pos + 2] = d100 % UInt8
        @inbounds buf[pos + 3] = (d100 >> 0x8) % UInt8
        pos += 4
    elseif pexp >= 100
        c = (pexp % 10) % UInt8
        @inbounds d100 = Base._dec_d100[div(pexp, 10) + 1]
        @inbounds buf[pos] = d100 % UInt8
        @inbounds buf[pos + 1] = (d100 >> 0x8) % UInt8
        @inbounds buf[pos + 2] = UInt8('0') + c
        pos += 3
    else
        @inbounds d100 = Base._dec_d100[pexp + 1]
        @inbounds buf[pos] = d100 % UInt8
        @inbounds buf[pos + 1] = (d100 >> 0x8) % UInt8
        pos += 2
    end
    return pos
end
