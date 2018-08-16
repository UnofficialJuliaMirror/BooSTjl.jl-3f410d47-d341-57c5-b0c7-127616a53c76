module BooSTjl

using StatsBase
using Base.LinAlg.BLAS

include("auxiliary.jl")
include("export.jl")

export Boost, BoostMore, estimate_derivatives, predictBoost

end # module
