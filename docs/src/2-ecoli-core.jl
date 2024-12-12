# # _E. coli_ core model

using ElementaryFluxModes

# It has been proven that enzyme constrained models with K constraints 
# will use a maximum of K EFMs in their optimal solution (de Groot 2019 (**check))

# Here, we take the _E. coli_ core model with two enzyme constraints, find the optimal 
# solution, and then find the EFMs in the optimal solution.
