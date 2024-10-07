"""
Package `ElementaryFluxModes` provides a Julia implementation of the Double Description 
method to calculate extreme rays of convex polyhedral cones. We follow the 
method described in Terzer 2009 thesis for the polyhedral cone Ρ = {x ∈ ℜ^d | Ax = 0, x >= 0}.
"""
module ElementaryFluxModes

using LinearAlgebra, RowEchelon
using DocStringExtensions

include("utils.jl")
export reorder_ns
export rational_nullspace
export make_all_irreversible

include("DDStandard.jl")
export DDStandard
export get_efms

include("DDBinary.jl")
export DDBinary
export make_bitmap
export check_adjacency

end
