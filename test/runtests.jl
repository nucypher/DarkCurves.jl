using Jute
using DarkCurves
using Random
using DarkIntegers

include("curve.test.jl")

exit(runtests(options=Dict(:exclude_tags => [:performance])))
