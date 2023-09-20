module ElementaryFluxModes

# Write your package code here.

using LinearAlgebra, RowEchelon

include("utils.jl") 
export make_all_irreversible
export reversible_EFMs
export fix_fluxes
export clean_DD_result
# export initialiseR
# export makeBitmap
export rational_nullspace
# export preprocessing
# export main_programme
# export cleanE
# export checkByLP

include("DDStandard.jl")
export DDStandard
export checkR
export adjacency_test
export DDStandardFixedFluxes

end
