# # Calculating and differentiating the optimal flux modes (OFMs) in a toy model

using ElementaryFluxModes
import COBREXA as X
import AbstractFBCModels as A
import Tulip as T
import FastDifferentiation as F
const Ex = F.Node
using DataFrames
# ## Load a simple model

# The code used to construct the model is located in `test/simple_model.jl`, but it is not shown here for brevity.

include("../../test/simple_model.jl"); #hide

model

N = A.stoichiometry(model)

# We wish to add an lower bound on the ATP maintenance proxy reaction, to illustrate how OFMs would work in a bigger model.

model.reactions["ATPM"].lower_bound = 10

rid_kcat
kcats = Symbol.(keys(rid_kcat))

float_reaction_isozymes = Dict{String,Dict{String,X.Isozyme}}()
for (rid, rxn) in model.reactions
    grrs = rxn.gene_association_dnf
    isnothing(grrs) && continue # skip if no grr available
    haskey(rid_kcat, rid) || continue # skip if no kcat data available
    for (i, grr) in enumerate(grrs)

        kcat = rid_kcat[rid] 

        d = get!(float_reaction_isozymes, rid, Dict{String,X.Isozyme}()) 
        d["isozyme_$i"] = X.Isozyme( 
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), 
            kcat_forward = kcat, 
            kcat_reverse = nothing,
        )
    end
end

gene_product_molar_masses

capacity

# ## Solve the toy model with enzyme constraints

ec_solution = X.enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes = float_reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer = T.Optimizer,
)

ec_solution.fluxes

@test all(values(ec_solution.fluxes) .> 1e-5) #src


# This optimal solution uses every reaction in the model, and since we have an enzyme capacity with two enzyme pools, theory has shown that there must be two OFMs in this optimal solution, and the optimal proportion of their flux is unique.

# To calculate these OFMs, we require the index and optimal value of ATPM and the objective reaction r6

atpm_idx = findfirst(x -> x == "ATPM", A.reactions(model))
biomass_idx = findfirst(x -> x == "r6", A.reactions(model))
fixed_fluxes = [atpm_idx, biomass_idx]
flux_values = [ec_solution.fluxes["ATPM"], ec_solution.fluxes["r6"]]

# Calculate OFMs

OFMs = get_ofms(Matrix(N), fixed_fluxes, flux_values)

@test OFMs ≈ [ #src
    10.0 10.0 #src
    75.0 0.0 #src
    75.0 75.0 #src
    75.0 0.0 #src
    75.0 0.0 #src
    0.0 75.0 #src
    65.0 65.0 #src
] || OFMs ≈ [ #src
    10.0 10.0 #src
    0.0 75.0 #src
    75.0 75.0 #src
    0.0 75.0 #src
    0.0 75.0 #src
    75 0.0 #src
    65.0 65.0 #src
] #src

OFM_dicts =
    [Dict(A.reactions(model) .=> OFMs[:, 1]), Dict(A.reactions(model) .=> OFMs[:, 2])]

# We see that the first OFM uses reactions 1,2,3,4,6 and ATPM, and the second OFM uses reactions 2,5,6 and ATPM.

OFM_dicts = [
    Dict(x => y / OFM_dicts[1]["r6"] for (x, y) in OFM_dicts[1]),
    Dict(x => y / OFM_dicts[2]["r6"] for (x, y) in OFM_dicts[2]),
]

# Calculate proportion of the OFMs used in the optimal solution

M = [
    OFM_dicts[1]["r5"] OFM_dicts[2]["r5"]
    OFM_dicts[1]["r1"] OFM_dicts[2]["r1"]
]

v = [
    ec_solution.fluxes["r5"]
    ec_solution.fluxes["r1"]
]

λ = M \ v

λ./sum(λ)

@test λ./sum(λ) ≈ [0.6666666667049761 0.33333333329502396]' #src
# Two thirds of the optimal flux is provided by the first OFM (using r3 and r4), and the remaining third by the second OFM (using r5).

# ## Differentiate the OFM usage

# It remains to differentiate to find the sensitivities of the OFM usage to the model parameters. For this, we must set up parameter isozymes.


parameter_values = Dict(
    Symbol(x) => iso.kcat_forward for (x, y) in float_reaction_isozymes for (_, iso) in y
)

parameters = Ex.(collect(keys(parameter_values)))
rid_gcounts = Dict(
    rid => [v.gene_product_stoichiometry for (k, v) in d][1] for
    (rid, d) in float_reaction_isozymes
)
rid_pid = Dict(rid => Ex(Symbol(rid)) for (rid, _) in rid_kcat)


sens = differentiate_ofm(
    OFM_dicts,
    parameters,
    rid_pid,
    parameter_values,
    rid_gcounts,
    capacity,
    gene_product_molar_masses,
    T.Optimizer,
)

# We may scale these sensitivities to calculate the control coefficients, (p/λ)*dλ/dp

control = Matrix(undef, size(sens, 1), size(sens, 2))
for (i, col) in enumerate(eachcol(sens))
    control[:, i] = collect(values(parameter_values))[i] .* col ./ λ
end
control

# We now use the order of `parameters` to put this data into a DataFrame, where each entry corresponds to the control coefficient (p/λ)*dλ/dp for p as the kcat of r3, r4 or r5, and for λ₁ and λ₂. We see that reaction 3 and 4 control the usage of the first OFM, whereas reaction 5 controls the usage of the second OFM. Increasing the kcat value of reaction 3 by 1% would increase the flux through OFM₁ by 0.5%, increasing the kcat of reaction 4 by 1% would have the same effect, and increasing the kcat of reaction 5 by 1% would increase the flux through OFM₂ by 1%.

df = DataFrame(
    Variable = ["λ₁", "λ₂"],
    r3 = control[:, 1],
    r4 = control[:, 3],
    r5 = control[:, 2],
)
