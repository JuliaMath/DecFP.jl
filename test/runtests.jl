using DecFP, Compat
using Base.Test

@test unsafe_load(DecFP.flags) == 0

import DecFP.isnanstr
@test isnanstr("nan") && isnanstr("  +NAN") && isnanstr("-NaN") && !isnanstr("nano")

for T in (Dec32, Dec64, Dec128)
    println(STDERR, "TESTING $T ...")

    if T == Dec32
        @test d32"3.2" * d32"4.5" == d32"14.4"
        @test eps(Dec32) == parse(Dec32, "1e-6")
        @test realmax(Dec32) == parse(Dec32, "+9999999E+90")
        @test realmin(Dec32) == parse(Dec32, "+1E-101")
        @test bswap(Dec32(1.5)) == reinterpret(Dec32, bswap(reinterpret(UInt32, Dec32(1.5))))
    elseif T == Dec64
        @test d"3.2" * d"4.5" == d"14.4"
        @test d64"3.2" * d64"4.5" == d64"14.4"
        @test eps(Dec64) == parse(Dec64, "1e-15")
        @test realmax(Dec64) == parse(Dec64, "+9999999999999999E+369")
        @test realmin(Dec64) == parse(Dec64, "+1E-398")
        @test bswap(Dec64(1.5)) == reinterpret(Dec64, bswap(reinterpret(UInt64, Dec64(1.5))))
    else
        @test d128"3.2" * d128"4.5" == d128"14.4"
        @test eps(Dec128) == parse(Dec128, "1e-33")
        @test realmax(Dec128) == parse(Dec128, "+9999999999999999999999999999999999E+6111")
        @test realmin(Dec128) == parse(Dec128, "+1E-6176")
        @test bswap(Dec128(1.5)) == reinterpret(Dec128, bswap(reinterpret(UInt128, Dec128(1.5))))
    end

    @test parse(T, "1.5")::T == 1.5
    @test isnan(parse(T, "NaN")::T) && isnan(T(NaN))
    @test isinf(parse(T, "inf")::T) && !isfinite(parse(T, "inf")::T)
    @test parse(T,"inf")::T == typemax(T)::T
    @test parse(T,"-inf")::T == typemin(T)::T

    @test parse(T, "0.1")::T == T(1//10)

    x,y,z = 1.5, -3.25, 0.0625 # exactly represented in binary
    xd = T(x); yd = T(y); zd = T(z)

    @test fma(xd,yd,zd)::T == muladd(xd,yd,zd)::T == xd*yd+zd == fma(x,y,z)

    @test one(T)::T == 1
    @test zero(T)::T == 0
    @test signbit(xd) == signbit(x) && signbit(yd) == signbit(y) && signbit(parse(T,"-Inf"))
    @test sign(xd)::T == sign(x) && sign(yd)::T == sign(y)

    @test nextfloat(xd)::T == xd + eps(T)
    @test prevfloat(xd)::T == xd - eps(T)
    @test nextfloat(T(1.5e10)) == 1.5e10 + eps(T(1.5e10))
    @test prevfloat(T(1.5e10)) == 1.5e10 - eps(T(1.5e10))

    for f in (isnan,isinf,isfinite,issubnormal,abs)
        @test f(xd) == f(x)
    end

    for f in (+,-,*,copysign,==,>,>=,<,<=)
        @test f(xd,yd) == f(x,y)
    end

    for f in (/,hypot,atan2,^)
        @test f(xd,yd) ≈ f(x,y)
        if f != ^
            @test f(yd,xd) ≈ f(y,x)
            @test f(yd,yd) ≈ f(y,y)
        end
    end

    for f in (exp,log,sin,cos,tan,asin,acos,atan,sinh,cosh,tanh,asinh,acosh,atanh,log1p,expm1,log10,log2,exp2,exp10,lgamma,sqrt,cbrt)
        v = f == acosh ? x : z
        @test f(T(v)) ≈ f(v)
    end

    for c in (π, e, γ, catalan, φ)
        @test T(c) ≈ Float64(c)
    end

    for Tf in (Float64, Float32, Float16)
        @test xd == Tf(x) == T(Tf(x)) == Tf(xd)
    end
    for Ti in (Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64)
        @test parse(T, "17") == T(Ti(17)) == Ti(17) == Ti(T(17))
        @test floor(Ti, T(4.5)) == 4 == ceil(Ti, T(4.5)) - 1
        @test_throws InexactError convert(Ti, xd)
        @test_throws InexactError floor(Ti, realmax(T))
        @test_throws InexactError ceil(Ti, realmax(T))
        if Ti <: Signed
            @test parse(T, "-17") == T(Ti(-17)) == Ti(-17) == Ti(T(-17))
            @test floor(Ti, T(-4.5)) == -5 == ceil(Ti, T(-4.5)) - 1
            @test_throws InexactError convert(Ti, yd)
        end
    end

    TI = eval(Symbol(string("UInt", sizeof(T)*8)))
    @test bswap(xd) == reinterpret(T, bswap(reinterpret(TI, xd)))

    @test_throws InexactError parse(T, "1e10000")
    @test_throws DomainError asin(xd)
    @test_throws DomainError sqrt(yd)
    @test_throws DomainError acosh(zd)

    # @test ldexp(parse(T, "1"), 3) == 1000
    # @test exponent(parse(T, "1000")) == 3
    @test sqrt(complex(yd)) ≈ sqrt(complex(y))

    @test typeof(xd * pi) == T
    @test typeof((xd+yd*im)*pi) == Complex{T}
end

@test unsafe_load(DecFP.flags) == 0
