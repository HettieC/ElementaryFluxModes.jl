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
