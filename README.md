# DecFP: IEEE Decimal Floating-point in Julia
[![Build Status](https://travis-ci.org/JuliaMath/DecFP.jl.svg)](https://travis-ci.org/JuliaMath/DecFP.jl) [![Build status](https://ci.appveyor.com/api/projects/status/si1d6og9wxsu8178?svg=true)](https://ci.appveyor.com/project/StevenGJohnson/decfp-jl) [![Coverage Status](https://coveralls.io/repos/github/JuliaMath/DecFP.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaMath/DecFP.jl?branch=master)

The DecFP package is a Julia wrapper around the [Intel Decimal
Floating-Point Math
Library](https://software.intel.com/en-us/articles/intel-decimal-floating-point-math-library),
providing a software implementation of the IEEE 754-2008 Decimal
Floating-Point Arithmetic specification.

32-bit, 64-bit, and 128-bit decimal floating-point types `Dec32`,
`Dec64`, and `Dec128`, respectively, are provided.  This is very
different from packages such as
[Decimals.jl](https://github.com/tinybike/Decimals.jl), which provide
*arbitrary-precision* decimal types analogous to `BigFloat`: arbitrary
precision types are very flexible, but fixed-precision types such
as those in DecFP are much faster (though still about 100x slower than
the hardware binary floating-point types `Float32` and `Float64`).

The DecFP package currently requires the Julia 0.5 release or later.

## Usage

`Dec64` and the other types mentioned above can be constructed from
other Julia numeric types (binary floating-point or integers) via
`Dec64(3.5)` or `Dec(3)`, from strings by `parse(Dec64, "3.2")` or
`d64"3.2"` (a Julia string macro); similarly for `Dec32` and `Dec128`.
The string macro `d"3.2"` constructs `Dec64`.

Once a decimal float is constructed, most Julia arithmetic and
special functions should work without modification.  For example,
`d"3.2" * d"4.5"` produces the `Dec64` result `+1440E-2` (14.4).
Most basic arithmetic functions are supported, and many special functions
(`sqrt`, `log`, trigonometric functions, etc.).   Mixed operations
involving decimal and binary floating-point or integer types are supported
(the result is promoted to decimal floating-point).

In general, you should be able to use the `DecFP` types in any context
where you would have used binary floating-point types: arrays, complex
arithmetic, and linear algebra should all work, for the most part.
