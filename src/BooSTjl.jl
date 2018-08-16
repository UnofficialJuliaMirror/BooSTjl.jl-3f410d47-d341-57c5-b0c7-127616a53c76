module BooSTjl

export Boost, BoostMore, estimate_derivatives, predictBoost

using StatsBase
using Base.LinAlg.BLAS


include("auxiliary.jl")
include("export.jl")

end # module
