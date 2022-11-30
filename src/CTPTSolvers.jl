module CTPTSolvers

using ForwardDiff
using LinearAlgebra
using StructArrays

export stability
export VanDerWaalsComponent
export VanDerWaalsMixture

include("types.jl")
include("constants.jl")
include("initials.jl")
include("nlsolve.jl")
include("solvecubic.jl")
include("stability.jl")

include("vanderwalls.jl")

end # module
