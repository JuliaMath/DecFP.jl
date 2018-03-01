using DecFP, Compat.Test, Compat.Printf

if !(VERSION < v"0.7.0-DEV.1592")
    using Base.MathConstants
end

@test unsafe_load(DecFP.flags[]) == 0

import DecFP.isnanstr
@test isnanstr("nan") && isnanstr("  +NAN") && isnanstr("-NaN") && !isnanstr("nano")

for T in (Dec32, Dec64, Dec128)
    info("TESTING $T ...")

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
    @test T("0.1")::T == T(1//10)

    io = IOBuffer()
    show(io, T("NaN")); @test String(take!(io)) == "NaN"
    show(io, T("Inf")); @test String(take!(io)) == "Inf"
    show(io, T("-Inf")); @test String(take!(io)) == "-Inf"
    show(io, T("0")); @test String(take!(io)) == "0.0"
    show(io, T("-0")); @test String(take!(io)) == "-0.0"
    show(io, T("1")); @test String(take!(io)) == "1.0"
    show(io, T("-1")); @test String(take!(io)) == "-1.0"
    show(io, T("1.000")); @test String(take!(io)) == "1.0"
    show(io, T("1e5")); @test String(take!(io)) == "100000.0"
    show(io, T("1e6")); @test String(take!(io)) == "1.0e6"
    show(io, T("1.23456e6")); @test String(take!(io)) == "1.23456e6"
    show(io, T("1e-1")); @test String(take!(io)) == "0.1"
    show(io, T("1e-4")); @test String(take!(io)) == "0.0001"
    show(io, T("1e-5")); @test String(take!(io)) == "1.0e-5"
    show(io, T("1.20e3")); @test String(take!(io)) == "1200.0"
    show(io, T("123.456")); @test String(take!(io)) == "123.456"
    show(io, T("0.00123456")); @test String(take!(io)) == "0.00123456"

    # some Dec128 tests fail due to Issue #47
    if T != Dec128
        @test @sprintf("%7.2f", T("1.2345")) == "   1.23"
        @test @sprintf("%-7.2f", T("1.2345")) == "1.23   "
        @test @sprintf("%07.2f", T("1.2345")) == "0001.23"
        @test @sprintf("%.0f", T("1.2345")) == "1"
        @test @sprintf("%#.0f", T("1.2345")) == "1."
        @test @sprintf("%.4e", T("1.2345")) == "1.2345e+00"
        @test @sprintf("%.4E", T("1.2345")) == "1.2345E+00"

        @test @sprintf("%6.2f", T("9.999")) == " 10.00"
        @test @sprintf("%9.2e", T("9.999e5")) == " 1.00e+06"

        @test @sprintf("%f", T("Inf")) == "Inf"
        @test @sprintf("%f", T("NaN")) == "NaN"

        @test @sprintf("%.0e", T("3e42")) == "3e+42"
        @test @sprintf("%#.0e", T("3e42")) == "3.e+42"

        @test @sprintf("%e", T("3e42")) == "3.000000e+42"
        @test @sprintf("%E", T("3e42")) == "3.000000E+42"
        @test @sprintf("%e", T("3e-42")) == "3.000000e-42"
        @test @sprintf("%E", T("3e-42")) == "3.000000E-42"

        @sprintf("%.6g", T("12345670.")) == "1.23457e+07"
        @sprintf("%.6g", T("1234567.")) == "1.23457e+06"
        @sprintf("%.6g", T("123456.7")) == "123457"
        @sprintf("%.6g", T("12345.67")) == "12345.7"
        @sprintf("%.6g", T("12340000.0")) == "1.234e+07"

        @test @sprintf("%10.5g", T("123.4")) == "     123.4"
        @test @sprintf("%+10.5g", T("123.4")) == "    +123.4"
        @test @sprintf("% 10.5g", T("123.4")) == "     123.4"
        @test @sprintf("%#10.5g", T("123.4")) == "    123.40"
        @test @sprintf("%-10.5g", T("123.4")) == "123.4     "
        @test @sprintf("%-+10.5g", T("123.4")) == "+123.4    "
        @test @sprintf("%010.5g", T("123.4")) == "00000123.4"
        @test @sprintf("%10.5g", T("-123.4")) == "    -123.4"
        @test @sprintf("%010.5g", T("-123.4")) == "-0000123.4"
        @test @sprintf("%.6g", T("12340000.0")) == "1.234e+07"
        @test @sprintf("%#.6g", T("12340000.0")) == "1.23400e+07"

        @test @sprintf("%.2f %.4f", T("12.34567"), 9.87) == "12.35 9.8700"
    end

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

    @test eps(maxintfloat(T) - 1) == 1
    @test eps(maxintfloat(T)) == 10
    @test maxintfloat(T) == maxintfloat(T)+1 > maxintfloat(T)-1 > maxintfloat(T)-2

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
    for Ti in (Integer,Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64)
        if Ti != Integer
            @test parse(T, "17") == T(Ti(17)) == Ti(17) == Ti(T(17))
        end
        @test trunc(Ti, T(2.7)) === floor(Ti, T(2.7)) === round(Ti, T(2.7), RoundDown) === round(Ti, T(2.7), RoundToZero) === Ti(2)
        @test ceil(Ti, T(2.3)) === round(Ti, T(2.3), RoundUp) === round(Ti, T(2.3), RoundFromZero) === Ti(3)
        @test round(Ti, T(1.5)) === round(Ti, T(2.5)) === round(Ti, T(1.5), RoundNearest) === round(Ti, T(2.5), RoundNearest) === Ti(2)
        @test round(Ti, T(2.5), RoundNearestTiesAway) === round(Ti, T(3.3), RoundNearestTiesAway) === Ti(3)
        @test round(Ti, T(2.5), RoundNearestTiesUp) === round(Ti, T(3.3), RoundNearestTiesUp) === Ti(3)
        @test_throws InexactError convert(Ti, xd)
        @test_throws InexactError trunc(Ti, realmax(T))
        @test_throws InexactError floor(Ti, realmax(T))
        @test_throws InexactError ceil(Ti, realmax(T))
        if Ti <: Signed
            @test parse(T, "-17") == T(Ti(-17)) == Ti(-17) == Ti(T(-17))
        end
        if Ti <: Signed || Ti === Integer
            @test trunc(Ti, T(-2.7)) === round(Ti, T(-2.7), RoundToZero) === Ti(-2)
            @test floor(Ti, T(-2.3)) === round(Ti, T(-2.3), RoundDown) === round(Ti, T(-2.3), RoundFromZero) === Ti(-3)
            @test ceil(Ti, T(-2.7)) === round(Ti, T(-2.7), RoundUp) === Ti(-2)
            @test round(Ti, T(-1.5)) === round(Ti, T(-2.5)) === round(Ti, T(-1.5), RoundNearest) === round(Ti, T(-2.5), RoundNearest) === Ti(-2)
            @test round(Ti, T(-2.5), RoundNearestTiesAway) === round(Ti, T(-3.3), RoundNearestTiesAway) === Ti(-3)
            @test round(Ti, T(-1.5), RoundNearestTiesUp) === round(Ti, T(-0.7), RoundNearestTiesUp) === Ti(-1)
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

@test unsafe_load(DecFP.flags[]) == 0

# issue #37
@test reinterpret(UInt128, Dec128(1.5)) == 0x303e000000000000000000000000000f
# issue #38
@test collect(v for i in 1:1, v in zeros(Dec128, 1)) ==  zeros(Dec128, 1, 1)
