# DecFP: IEEE Decimal Floating-point in Julia

[![CI](https://github.com/JuliaMath/DecFP.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaMath/DecFP.jl/actions/workflows/CI.yml)
[![Coverage Status](https://coveralls.io/repos/github/JuliaMath/DecFP.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaMath/DecFP.jl?branch=master)
[![codecov](https://codecov.io/github/JuliaMath/DecFP.jl/graph/badge.svg?token=9thK9TVT9V)](https://codecov.io/github/JuliaMath/DecFP.jl)

The DecFP package is a Julia wrapper around the [Intel Decimal
Floating-Point Math
Library](https://software.intel.com/en-us/articles/intel-decimal-floating-point-math-library),
providing a software implementation of the IEEE 754-2008 Decimal
Floating-Point Arithmetic specification.

32-bit, 64-bit, and 128-bit decimal floating-point types `Dec32`,
`Dec64`, and `Dec128`, respectively, are provided (corresponding to 7, 16, and 34 decimal [digits of precision](https://en.wikipedia.org/wiki/Significand), respectively).  This is very
different from packages such as
[Decimals.jl](https://github.com/tinybike/Decimals.jl), which provide
*arbitrary-precision* decimal types analogous to `BigFloat`: arbitrary
precision types are very flexible, but fixed-precision types such
as those in DecFP are much faster (though still about 50x slower than
the hardware binary floating-point types `Float32` and `Float64`) and
more memory-efficient (an array of `Dec64` values has exactly the
same memory footprint as an array of `Float64` values).

The latest version of the DecFP package requires Julia 1.7 or later.

## Usage

`Dec64` and the other types mentioned above can be constructed from
other Julia numeric types (binary floating-point or integers) via
`Dec64(3.5)` or `Dec(3)`, from strings by `parse(Dec64, "3.2")` or
`d64"3.2"` (a Julia string macro); similarly for `Dec32` and `Dec128`.
The string macro `d"3.2"` constructs `Dec64`.

* **Note:** A decimal floating-point constant like `1.1` in Julia refers to a `Float64` (binary floating-point) value.  So, for example, `Dec128(1.999999999999999999) == 2.0 ≠ d128"1.999999999999999999"`, since Julia first rounds `1.999999999999999999` to the closest `Float64` value (`2.0`) before converting to `Dec128`.  If you need to specify an exact decimal constant, therefore, use one of the string-based constructors.  If you use a string macro like `d128"1.999999999999999999"`, then the string parsing occurs *before* compilation and incurs no runtime cost.

Once a decimal float is constructed, most Julia arithmetic and
special functions should work without modification.  For example,
`d"3.2" * d"4.5"` produces the (exact) `Dec64` result `14.4`
All basic arithmetic functions are supported, and many special functions
(`sqrt`, `log`, trigonometric functions, etc.).   Mixed operations
involving decimal and binary floating-point or integer types are supported
(the result is promoted to decimal floating-point).

In general, you should be able to use the `DecFP` types in any context
where you would have used binary floating-point types: arrays, complex
arithmetic, and linear algebra should all work, for the most part.
