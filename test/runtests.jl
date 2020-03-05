using Jute
using Random
using DarkIntegers
using BenchmarkTools

using DarkCurves

include("TestUtils.jl")
include("curve.test.jl")
include("endomorphism.test.jl")

exit(runtests(options=Dict(:exclude_tags => [:performance])))
