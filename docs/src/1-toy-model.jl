# # Calculating EFMs in a toy model

using ElementaryFluxModes

import AbstractFBCModels as A

# ## Load a simple model

# The code used to construct the model is located in `test/simple_model.jl`, but
# it is not shown here for brevity.

include("../../test/simple_model.jl"); #hide

model

model.reactions

delete!(model.reactions, "ATPM") #src

# Get the stoichiometric matrix of the model, this is what we use to find EFMs

N = A.stoichiometry(model)

# ## Calculate the EFMs

# Run the double description algorith on `N` and `K`

E = get_efms(Matrix(N))

# If preferred, we can transform the vector of EFMs, **`E`**, into a dictionary
# of reaction => fluxes through the efms

EFM_dict = Dict(A.reactions(model) .=> eachrow(E))
