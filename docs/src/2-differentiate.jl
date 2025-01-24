# # Differentiating the EFMs of a toy model

# The optimal flux distribution of any metabolic model can be written as a 
# weighted sum of the EFMs of that model. We are interested in calculating 
# the sensitivity of these weightings to the model parameters, and can use 
# DifferentiableMetabolism.jl to efficiently do so.

using ElementaryFluxModes

import AbstractFBCModels as A
import FastDifferentiation as F
const Ex = F.Node
import DifferentiableMetabolism as D
using COBREXA
import COBREXA as X
using JSON
import Tulip as T

# ## Build a simple enzyme constrained model

# The code used to construct the model is located in `test/simple_model.jl`, but
# it is not shown here for brevity.

include("../../test/simple_model.jl"); #hide

model

parameter_values = Dict{Symbol,Float64}()
reaction_isozymes = Dict{String,Dict{String,X.IsozymeT{Ex}}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
float_reaction_isozymes = Dict{String,Dict{String,X.Isozyme}}() #src
for rid in A.reactions(model)
    grrs = A.reaction_gene_association_dnf(model, rid)
    isnothing(grrs) && continue # skip if no grr available
    haskey(rid_kcat, rid) || continue # skip if no kcat data available
    for (i, grr) in enumerate(grrs)

        kcat = rid_kcat[rid]
        parameter_values[Symbol(rid)] = kcat

        d = get!(reaction_isozymes, rid, Dict{String,X.IsozymeT{Ex}}())
        d["isozyme_$i"] = X.IsozymeT{Ex}(
            gene_product_stoichiometry = Dict(grr .=> fill(float(1.0), size(grr))), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = Ex(Symbol(rid)),
            kcat_reverse = 0.0,
        )
        d2 = get!(float_reaction_isozymes, rid, Dict{String,X.Isozyme}()) #src
        d2["isozyme_$i"] = X.Isozyme( #src
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), #src
            kcat_forward = kcat, #src
            kcat_reverse = 0.0, #src
        ) #src
    end
end

km = enzyme_constrained_flux_balance_constraints( # kinetic model
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
)

ec_solution = D.optimized_values(
    km,
    parameter_values;
    objective = km.objective.value,
    optimizer = T.Optimizer,
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

ec_solution.tree.fluxes

ec_solution_fba = enzyme_constrained_flux_balance_analysis( #src
    model; #src
    reaction_isozymes = float_reaction_isozymes, #src
    gene_product_molar_masses, #src
    capacity, #src
    optimizer = T.Optimizer, #src
) #src

@test isapprox(ec_solution.objective, ec_solution_fba.objective; atol = TEST_TOLERANCE) #src

@test any(isapprox.(values(ec_solution.gene_product_amounts), 0, atol = 1e-8)) #src

# We have a solution that uses every reaction, and the enzyme capacities are both full.
# Therefore, we may calculate the EFMs of this solution and directly differentiate
# them, with no pruning required.

# ## Calculate EFMs of the optimal solution

# We need to input the stoichiometric matrix `N` into ElementaryFluxModes.jl

N = A.stoichiometry(model)

# Calculate a flux matrix of the EFMs, the size of which is (n,k), for n reactions 
# and k EFMs

E = get_efms(Matrix(N))

# Make a dictionary out of the EFM result

EFM_dict = Dict(A.reactions(model) .=> eachrow(E))
EFMs = [
    Dict(k => v[1] / EFM_dict["r6"][1] for (k, v) in EFM_dict),
    Dict(k => v[2] / EFM_dict["r6"][2] for (k, v) in EFM_dict),
]

@test EFM_dict == Dict(
    "r1" => [1.0, 0.0],
    "r2" => [1.0, 1.0],
    "r5" => [0.0, 1.0],
    "r6" => [1.0, 1.0],
    "r3" => [1.0, 0.0],
    "r4" => [1.0, 0.0]) #src

# The optimal solution, **v**, can be written as λ₁**EFM₁**+λ₂**EFM₂**=**v**
# so that the λ give us the weightings of the two EFMs.

# Let's calculate λ₁ and λ₂, using reactions `r1` and `r5`, as these are not 
# shared by the EFMs

M = [
    EFM_dict["r1"][1] EFM_dict["r1"][2]
    EFM_dict["r5"][1] EFM_dict["r5"][2]
]

v = [
    ec_solution.tree.fluxes["r1"]
    ec_solution.tree.fluxes["r5"]
]

λ = M \ v

# The optimal solution is therefore made up of 50 units of flux through EFM₁
# and 25 units of flux through EFM₂.

# ## Differentiate the EFMs

# We have calculated the EFMs, and now wish to differentiate their weightings, 
# `λ`, with respect to the model parameters

parameters = Ex.(collect(keys(parameter_values)))
p_vals = collect(values(parameter_values))
rid_pid =
    Dict(rid => [iso.kcat_forward for (k, iso) in v][1] for (rid, v) in reaction_isozymes)

sens_efm = differentiate_efm(
    EFMs,
    parameters,
    rid_pid,
    parameter_values,
    rid_gcounts,
    capacity,
    gene_product_molar_masses,
    T.Optimizer,
)

@test sens_efm == [
    2.5 -0.0 2.5
    -0.0 5.0 -0.0
] #src 
