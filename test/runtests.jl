using DecFP, Test, Printf, Base.MathConstants, SpecialFunctions

@test DecFP.flags[Threads.threadid()] == 0

import DecFP.isnanstr
@test isnanstr("nan") && isnanstr("  +NAN") && isnanstr("-NaN") && !isnanstr("nano")

function testthreads(T, i, mode)
    @test @sprintf("%.0f", T(i)) == string(i)
    @test rounding(T) == RoundNearest
    setrounding(T, mode) do
        @test rounding(T) == mode
        if mode in (RoundNearest, RoundToZero, RoundDown)
            @test T(1) + eps(T(1)) / 2 == T(1)
        else
            @test T(1) + eps(T(1)) / 2 == T(1) + eps(T(1))
        end
    end
    @test rounding(T) == RoundNearest
    @test_throws InexactError convert(Int64, T("1.5"))
    @test_throws DomainError sqrt(T(-2))
end

for T in (Dec32, Dec64, Dec128)
    @info "TESTING $T    nthreads = $(Threads.nthreads()) ..."

    if T == Dec32
        @test d32"3.2" * d32"4.5" == d32"14.4"
        @test T(1, Int32(maxintfloat(T) - 1), 90) == floatmax(T)
        @test_throws InexactError T(1, Int32(maxintfloat(T)), 90)
        @test_throws InexactError T(1, Int32(maxintfloat(T)) + 1, 0)
        @test T(1, 1, 96) == T("1e96")
        @test T(1, 1, -101) == nextfloat(T(0))
        @test T(1, 10, -102) == nextfloat(T(0))
        @test_throws InexactError T(1, 1, -102)
        @test_throws InexactError T(1, 11, -102)
        @test eps(Dec32) == parse(Dec32, "1e-6")
        @test floatmax(Dec32) == parse(Dec32, "+9999999E+90")
        @test bswap(Dec32(1.5)) == reinterpret(Dec32, bswap(reinterpret(UInt32, Dec32(1.5))))
    elseif T == Dec64
        @test d"3.2" * d"4.5" == d"14.4"
        @test d64"3.2" * d64"4.5" == d64"14.4"
        @test T(1, Int64(maxintfloat(T) - 1), 369) == floatmax(T)
        @test_throws InexactError T(1, Int64(maxintfloat(T)), 369)
        @test_throws InexactError T(1, Int64(maxintfloat(T)) + 1, 0)
        @test T(1, 1, 384) == T("1e384")
        @test T(1, 1, -398) == nextfloat(T(0))
        @test T(1, 10, -399) == nextfloat(T(0))
        @test_throws InexactError T(1, 1, -399)
        @test_throws InexactError T(1, 11, -399)
        @test eps(Dec64) == parse(Dec64, "1e-15")
        @test floatmax(Dec64) == parse(Dec64, "+9999999999999999E+369")
        @test bswap(Dec64(1.5)) == reinterpret(Dec64, bswap(reinterpret(UInt64, Dec64(1.5))))
    else
        @test d128"3.2" * d128"4.5" == d128"14.4"
        @test T(1, Int128(maxintfloat(T)) - 1, 6111) == floatmax(T)
        @test_throws InexactError T(1, Int128(maxintfloat(T)), 6111)
        @test_throws InexactError T(1, Int128(maxintfloat(T)) + 1, 0)
        @test T(1, 1, 6144) == T("1e6144")
        @test T(1, 1, -6176) == nextfloat(T(0))
        @test T(1, 10, -6177) == nextfloat(T(0))
        @test_throws InexactError T(1, 1, -6177)
        @test_throws InexactError T(1, 11, -6177)
        @test eps(Dec128) == parse(Dec128, "1e-33")
        @test floatmax(Dec128) == parse(Dec128, "+9999999999999999999999999999999999E+6111")
        @test bswap(Dec128(1.5)) == reinterpret(Dec128, bswap(reinterpret(UInt128, Dec128(1.5))))
    end

    @test parse(T, "1.5")::T == 1.5
    @test isnan(parse(T, "NaN")::T) && isnan(T(NaN))
    @test isinf(parse(T, "inf")::T) && !isfinite(parse(T, "inf")::T)
    @test parse(T,"inf")::T == typemax(T)::T
    @test parse(T,"-inf")::T == typemin(T)::T
    @test_throws ArgumentError parse(T, "1.0.1")
    @test_throws ArgumentError parse(T, "NaNny")
    @test tryparse(T, "1.5")::T == 1.5
    @test isnan(tryparse(T, "NaN"))
    @test tryparse(T, "NaNny") === nothing
    @test tryparse(T, "1.0.1") === nothing

    @test parse(T, "0.1")::T == T(1//10)
    @test T("0.1")::T == T(1//10)

    if T != Dec128
        @test T(0.1) == T("0.1")
        @test T(0.1, RoundNearest) == T("0.1")
        @test T(0.1, RoundDown) == T("0.1")
        @test T(0.1, RoundUp) == nextfloat(T("0.1"))
    end
    @test T("1", RoundDown) == T(1)
    @test T("1", RoundUp) == T(1)
    @test T("1.0000000000000000000000000000000000000001") == T(1)
    @test T("1.0000000000000000000000000000000000000001", RoundNearest) == T(1)
    @test T("1.0000000000000000000000000000000000000001", RoundDown) == T(1)
    @test T("1.0000000000000000000000000000000000000001", RoundUp) == nextfloat(T(1))
    @test T("0.9999999999999999999999999999999999999999") == T(1)
    @test T("0.9999999999999999999999999999999999999999", RoundNearest) == T(1)
    @test T("0.9999999999999999999999999999999999999999", RoundUp) == T(1)
    @test T("0.9999999999999999999999999999999999999999", RoundDown) == prevfloat(T(1))

    @test T(0, 0, 0) == T(0)
    @test T(1, 0, 10) == T(0)
    @test isequal(T(-1, 0, 0), T(-0.0))
    @test T(1, 125, -2) == T(1.25)
    @test T(-1, 500, -2) == T(-5)
    @test T(1, 5, 2) == T(500)
    @test_throws DomainError T(0, 1, 0)
    @test_throws DomainError T(2, 1, 0)
    @test T(125, -2) == T(1.25)
    @test T(-125, -2) == T(-1.25)
    @test T(1, Int128(maxintfloat(T)), 0) == maxintfloat(T)
    @test T(1, Int128(maxintfloat(T)) * 100, -2) == maxintfloat(T)

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

    x,y,z = 1.5, -3.25, 0.0625 # exactly represented in binary
    xd = T(x); yd = T(y); zd = T(z)

    @test fma(xd,yd,zd)::T == muladd(xd,yd,zd)::T == xd*yd+zd == fma(x,y,z)

    @test one(T)::T == 1
    @test zero(T)::T == 0
    @test signbit(xd) == signbit(x) && signbit(yd) == signbit(y) && signbit(parse(T,"-Inf"))
    @test isnan(sign(T(NaN)))
    @test isequal(sign(T(0.0)), T(0.0))
    @test isequal(sign(T(-0.0)), T(-0.0))
    @test sign(xd)::T == sign(x) && sign(yd)::T == sign(y)

    @test nextfloat(xd)::T == xd + eps(T)
    @test prevfloat(xd)::T == xd - eps(T)
    @test nextfloat(T(1.5e10)) == 1.5e10 + eps(T(1.5e10))
    @test prevfloat(T(1.5e10)) == 1.5e10 - eps(T(1.5e10))

    @test floatmin(T) > 0
    @test !issubnormal(floatmin(T))
    @test issubnormal(prevfloat(floatmin(T)))
    @test floatmin(T) == nextfloat(zero(T)) / eps(T)

    @test eps(maxintfloat(T) - 1) == 1
    @test eps(maxintfloat(T)) == 10
    @test maxintfloat(T) == maxintfloat(T)+1 > maxintfloat(T)-1 > maxintfloat(T)-2

    for f in (isnan,isinf,isfinite,issubnormal,abs)
        @test f(xd) == f(x)
    end

    for f in (+,-,*,copysign,==,>,>=,<,<=)
        @test f(xd,yd) == f(x,y)
    end

    for f in (/,hypot,atan,^)
        @test f(xd,yd) ≈ f(x,y)
        if f != ^
            @test f(yd,xd) ≈ f(y,x)
            @test f(yd,yd) ≈ f(y,y)
        end
    end

    for f in (exp,log,sin,cos,tan,asin,acos,atan,sinh,cosh,tanh,asinh,acosh,atanh,log1p,expm1,log10,log2,exp2,exp10,sqrt,cbrt)
        v = f == acosh ? x : z
        @test f(T(v)) ≈ f(v)
    end

    # issue #47
    @test exp10(T(0)) == T(1)

    for x = -3.0:0.25:3.0
        @test all(logabsgamma(T(x)) .≈ logabsgamma(x))
    end
    for x in (NaN, -Inf, -0.0, Inf)
        @test isequal(logabsgamma(T(x)), logabsgamma(x))
    end

    @test_throws DomainError gamma(T(-1))
    @test isnan(gamma(T(NaN)))
    for x in (-1.5, -0.5, -0.0, 0.0, 0.5, 1.0, 1.5, 2.0, Inf)
        @test gamma(T(x)) ≈ gamma(x)
    end

    for c in (π, e, γ, catalan, φ)
        @test T(c) ≈ Float64(c)
    end

    for Tf in (Float64, Float32, Float16)
        @test xd == Tf(x) == T(Tf(x)) == Tf(xd)
    end

    # issue #92
    p = -2
    @test T(2)^-2 == parse(T, "0.25") == T(2)^p

    # exercise literal_pow
    @test T(2)^0 === T(1)
    @test T(2)^1 === T(2)
    @test T(2)^2 === T(4)
    @test T(2)^3 === T(8)

    @test T("0.1") == 1//10
    @test 1//10 == T("0.1")
    @test T("-0.1") == -1//10
    @test T("0.2") == 1//5
    @test T("0.5") == 1//2
    @test T(1) == 1//1
    @test T(0) == 0//1
    @test T(Inf) == 1//0
    @test T(-Inf) == -1//0
    @test T(1) / T(3) != 1//3
    @test T(1) / T(3) < 1//3
    @test T(1) / T(3) <= 1//3
    @test T(2) / T(3) > 2//3
    @test T(2) / T(3) >= 2//3
    @test T(1) / T(6) != 1//6
    @test T(7) / T(2) == 7//2
    @test T(7) / T(5) == 7//5
    @test T(7) / T(100) == 7//100
    @test T(7) / T(300) != 7//300

    for Tf in (Float64, Float32, Float16)
        @test parse(T, "0.1") != parse(Tf, "0.1")
        @test parse(Tf, "0.1") != parse(T, "0.1")
        @test parse(T, "0.7") != parse(Tf, "0.7")
        @test parse(Tf, "0.7") != parse(T, "0.7")
        if Tf == Float16
            @test parse(T, "0.1") > parse(Tf, "0.1")
            @test parse(T, "0.1") >= parse(Tf, "0.1")
            @test parse(Tf, "0.1") < parse(T, "0.1")
            @test parse(Tf, "0.1") <= parse(T, "0.1")
            @test parse(T, "0.7") < parse(Tf, "0.7")
            @test parse(T, "0.7") <= parse(Tf, "0.7")
            @test parse(Tf, "0.7") > parse(T, "0.7")
            @test parse(Tf, "0.7") >= parse(T, "0.7")
        else
            @test parse(T, "0.1") < parse(Tf, "0.1")
            @test parse(T, "0.1") <= parse(Tf, "0.1")
            @test parse(Tf, "0.1") > parse(T, "0.1")
            @test parse(Tf, "0.1") >= parse(T, "0.1")
            @test parse(T, "0.7") > parse(Tf, "0.7")
            @test parse(T, "0.7") >= parse(Tf, "0.7")
            @test parse(Tf, "0.7") < parse(T, "0.7")
            @test parse(Tf, "0.7") <= parse(T, "0.7")
        end
        @test parse(T, "0.5") == parse(Tf, "0.5")
        @test parse(T, "0.5") <= parse(Tf, "0.5")
        @test parse(T, "0.5") >= parse(Tf, "0.5")
        @test parse(Tf, "0.5") == parse(T, "0.5")
        @test parse(Tf, "0.5") <= parse(T, "0.5")
        @test parse(Tf, "0.5") >= parse(T, "0.5")
    end

    @test hash(T(NaN)) == hash(NaN)
    @test hash(T(Inf)) == hash(Inf)
    @test hash(T(-Inf)) == hash(-Inf)
    @test hash(T(0)) == hash(0.0) == hash(0)
    @test hash(T(-0.0)) == hash(-0.0)
    @test hash(T(1)) == hash(1.0)
    @test hash(T(-1)) == hash(-1.0)
    @test hash(T(25)) == hash(25.0)
    @test hash(T(-25)) == hash(-25.0)
    @test hash(T(77)) == hash(77.0)
    @test hash(T(-77)) == hash(-77.0)
    @test hash(T(3500)) == hash(3500.0)
    @test hash(T(-3500)) == hash(-3500.0)
    @test hash(T(1234.5)) == hash(1234.5)
    @test hash(T(-1234.5)) == hash(-1234.5)
    @test hash(T(123.0625)) == hash(123.0625)
    @test hash(T(-123.0625)) == hash(-123.0625)
    @test hash(T(0.0009765625)) == hash(0.0009765625)
    @test hash(T(-0.0009765625)) == hash(-0.0009765625)
    @test hash(T(10)^70) == hash(BigInt(10)^70)

    for x = -5.0:0.25:5.0, y = -5.0:0.25:5.0
        @test isequal(rem(T(x), T(y)), rem(x, y))
        @test isequal(rem(T(x), T(y), RoundNearest), rem(x, y, RoundNearest))
        @test isequal(rem(T(x), T(y), RoundToZero), rem(x, y, RoundToZero))
        @test isequal(rem(T(x), T(y), RoundDown), rem(x, y, RoundDown))
        @test isequal(rem(T(x), T(y), RoundUp), rem(x, y, RoundUp))
        @test isequal(mod(T(x), T(y)), mod(x, y))
        @test isequal(T(x) % T(y), x % y)
    end
    for x in (NaN, -Inf, -0.0, 0.0, 1.0, Inf), y in (NaN, -Inf, -0.0, 0.0, 1.0, Inf)
        @test isequal(rem(T(x), T(y)), rem(x, y))
        @test isequal(rem(T(x), T(y), RoundNearest), rem(x, y, RoundNearest))
        @test isequal(rem(T(x), T(y), RoundToZero), rem(x, y, RoundToZero))
        @test isequal(rem(T(x), T(y), RoundDown), rem(x, y, RoundDown))
        @test isequal(rem(T(x), T(y), RoundUp), rem(x, y, RoundUp))
        @test isequal(mod(T(x), T(y)), mod(x, y))
        @test isequal(T(x) % T(y), x % y)
    end
    for x = -5.0:0.25:5.0
        @test isequal(modf(T(x)), modf(x))
    end
    for x in (NaN, -Inf, -0.0, Inf)
        @test isequal(modf(T(x)), modf(x))
    end

    @test trunc(T(2.7)) === floor(T(2.7)) === round(T(2.7), RoundDown) === round(T(2.7), RoundToZero) === T(2)
    @test ceil(T(2.3)) === round(T(2.3), RoundUp) === round(T(2.3), RoundFromZero) === T(3)
    @test round(T(1.5)) === round(T(2.5)) === round(T(1.5), RoundNearest) === round(T(2.5), RoundNearest) === T(2)
    @test round(T(2.5), RoundNearestTiesAway) === round(T(3.3), RoundNearestTiesAway) === T(3)
    @test round(T(2.5), RoundNearestTiesUp) === round(T(3.3), RoundNearestTiesUp) === T(3)

    for Ti in (Integer,Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64,Int128,UInt128)
        if Ti != Integer
            @test parse(T, "17") == T(Ti(17)) == Ti(17) == Ti(T(17))
        end
        @test trunc(Ti, T(2.7)) === floor(Ti, T(2.7)) === round(Ti, T(2.7), RoundDown) === round(Ti, T(2.7), RoundToZero) === Ti(2)
        @test ceil(Ti, T(2.3)) === round(Ti, T(2.3), RoundUp) === round(Ti, T(2.3), RoundFromZero) === Ti(3)
        @test round(Ti, T(1.5)) === round(Ti, T(2.5)) === round(Ti, T(1.5), RoundNearest) === round(Ti, T(2.5), RoundNearest) === Ti(2)
        @test round(Ti, T(2.5), RoundNearestTiesAway) === round(Ti, T(3.3), RoundNearestTiesAway) === Ti(3)
        @test round(Ti, T(2.5), RoundNearestTiesUp) === round(Ti, T(3.3), RoundNearestTiesUp) === Ti(3)
        @test_throws InexactError convert(Ti, xd)
        @test_throws InexactError trunc(Ti, floatmax(T))
        @test_throws InexactError floor(Ti, floatmax(T))
        @test_throws InexactError ceil(Ti, floatmax(T))
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

    @test rounding(T) == RoundNearest
    for mode in (RoundNearest, RoundToZero, RoundFromZero, RoundUp, RoundDown)
        setrounding(T, mode) do
            @test rounding(T) == mode
        end
    end
    @test rounding(T) == RoundNearest

    TI = eval(Symbol(string("UInt", sizeof(T)*8)))
    @test bswap(xd) == reinterpret(T, bswap(reinterpret(TI, xd)))

    @test parse(T, "1e10000") == T(Inf)
    @test_throws DomainError asin(xd)
    @test_throws DomainError sqrt(yd)
    @test_throws DomainError acosh(zd)

    @test ldexp10(parse(T, "1"), 3) == 1000
    @test exponent10(parse(T, "1000")) == 3
    @test sqrt(complex(yd)) ≈ sqrt(complex(y))

    @test typeof(xd * pi) == T
    @test typeof((xd+yd*im)*pi) == Complex{T}

    Threads.@threads for (i, mode) in collect(enumerate((RoundNearest, RoundToZero, RoundFromZero, RoundUp, RoundDown)))
        testthreads(T, i, mode)
    end

    # issue #85
    @test T(1.5) == T(T(1.5))
end

@test DecFP.flags[Threads.threadid()] == 0

@test widen(Dec32) == Dec64
@test widen(Dec64) == Dec128
@test widen(Dec32(1)) === Dec64(1)
@test widen(Dec64(1)) === Dec128(1)

@test hash(Dec32(9999999)) == hash(9999999)
@test hash(Dec32(-9999999)) == hash(-9999999)
@test hash(Dec64(9999999999999999)) == hash(9999999999999999)
@test hash(Dec64(-9999999999999999)) == hash(-9999999999999999)
@test hash(Dec128("9999999999999999999999999999999999")) == hash(9999999999999999999999999999999999)
@test hash(Dec128("-9999999999999999999999999999999999")) == hash(-9999999999999999999999999999999999)
@test hash(Dec128(2)^-30) == hash(BigFloat(2)^-30)
@test hash(-Dec128(2)^-30) == hash(-BigFloat(2)^-30)

# issue #37
@test reinterpret(UInt128, Dec128(1.5)) == 0x303e000000000000000000000000000f
# issue #38
@test collect(v for i in 1:1, v in zeros(Dec128, 1)) ==  zeros(Dec128, 1, 1)

@test setrounding(Dec64, RoundDown) do
    Float64(d64"1e100") < 1e100
end

@test Float64(d64"1e100") == 1e100

# issue #79
@test d64"1e100" != 1e100

# issue #93
@test parse(Dec64, "3935767060.093896713") == d64"3.935767060093897e9" ==
      Dec64(d128"3935767060.093896713")
