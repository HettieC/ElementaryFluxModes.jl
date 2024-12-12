# # Toy Model

using ElementaryFluxModes

using LinearAlgebra

# ## Load a simple model

# The code used to construct the model is located in `test/simple_model.jl`, but
# it is not shown here for brevity.

include("../../test/simple_model.jl"); #hide


model

# Get the stoichiometric matrix of the model, this is what we use to find EFMs

N = AbstractFBCModels.stoichiometry(model)

# ## Calculate the EFMs 

# Run the double description algorith on `N` and `K`

E = get_efms(Matrix(N))

# If preferred, we can transform the vector of EFMs, **`E`**, into a dictionary 
# of reaction => fluxes through the efms

EFM_dict = Dict(AbstractFBCModels.reactions(model) .=> eachrow(E))

