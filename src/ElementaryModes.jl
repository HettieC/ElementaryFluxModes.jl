module ElementaryModes

# Write your package code here.

using LinearAlgebra, RowEchelon, JuMP

include("utils.jl")
export initialiseR
export makeBitmap
export rational_nullspace
export preprocessing
export main_programme
export cleanE
export checkByLP

include("DDMethod.jl")
export DDStandard
export checkR

end
