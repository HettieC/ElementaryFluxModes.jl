# # _E. coli_ core model

using ElementaryFluxModes

import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
using JSONFBCModels
import FastDifferentiation as F
const Ex = F.Node
import DifferentiableMetabolism as D
using COBREXA
import COBREXA as X
using JSON
import Tulip as T

# It has been proven that enzyme constrained models with K constraints 
# will use a maximum of K EFMs in their optimal solution (de Groot 2019).

# Here, we take the _E. coli_ core model with two enzyme constraints, find the optimal 
# solution, and then find the EFMs in the optimal solution.


!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

include("../../test/data_static.jl")

model = convert(CM.Model, X.load_model("e_coli_core.json"))
# unconstrain glucose
model.reactions["EX_glc__D_e"].lower_bound = -1000.0
# constrain PFL to zero
model.reactions["PFL"].upper_bound = 0.0

rid_kcat = Dict(k => Ex(Symbol(k)) for (k, _) in ecoli_core_reaction_kcats)
kcats = Symbol.(keys(ecoli_core_reaction_kcats))

float_reaction_isozymes = Dict{String,Dict{String,X.Isozyme}}()
for (rid, rxn) in model.reactions
    grrs = rxn.gene_association_dnf
    isnothing(grrs) && continue # skip if no grr available
    haskey(ecoli_core_reaction_kcats, rid) || continue # skip if no kcat data available
    for (i, grr) in enumerate(grrs)

        kcat = ecoli_core_reaction_kcats[rid] * 3.6 # change unit to k/h

        d = get!(float_reaction_isozymes, rid, Dict{String,X.Isozyme}())
        d["isozyme_$i"] = X.Isozyme(
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = kcat, # assume forward and reverse have the same kcat
            kcat_reverse = kcat,
        )
    end
end

gene_product_molar_masses = Dict(k => v for (k, v) in ecoli_core_gene_product_masses)

capacity = [
    (
        "membrane", # enzyme group names must be of type Symbol
        [x for (x, y) in ecoli_core_subcellular_location if y == "membrane"],
        15.0,
    ),
    ("cytosol", [x for (x, y) in ecoli_core_subcellular_location if y == "cytosol"], 35.0),
]

ec_solution = X.enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes = float_reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer = T.Optimizer,
)

ec_solution

# This solution is both producing acetate and consuming oxygen, therefore it looks like 
# overflow metabolism.

ec_solution.fluxes["EX_ac_e"]
ec_solution.fluxes["EX_o2_e"]

# Since we had two enzyme pools, the membrane and the cytosol, this optimal 
# solution will be composed of at most two EFMs. In order to find out how 
# small parameter changes will affect this optimal composition, we must first 
# remove inactive reactions, so that the optimal solution is a unique 
# superposition of the two EFMs, and we can differentiate it.

# This solution contains many inactive reactions
sort(collect(ec_solution.fluxes), by = ComposedFunction(abs, last))

@test any(values(ec_solution.fluxes) .≈ 0) #src

# And also many inactive gene products. 

sort(collect(ec_solution.gene_product_amounts), by = last)

@test any(isapprox.(values(ec_solution.gene_product_amounts), 0, atol = 1e-8)) #src

# Let us start by pruning the model.

flux_zero_tol = 1e-6 # these bounds make a real difference!
gene_zero_tol = 1e-6
pruned_model, pruned_reaction_isozymes = D.prune_model(
    model,
    ec_solution.fluxes,
    ec_solution.gene_product_amounts,
    float_reaction_isozymes,
    ec_solution.isozyme_forward_amounts,
    ec_solution.isozyme_reverse_amounts,
    flux_zero_tol,
    gene_zero_tol,
)

# Since we will later be differentiating our solution, we need to make 
# parameter isozymes and a kinetic model

parameter_values = Dict(
    Symbol(x) => iso.kcat_forward for (x, y) in pruned_reaction_isozymes for (_, iso) in y
)

rid_kcat = Dict(k => Ex(Symbol(k)) for (k, _) in parameter_values)
parameter_isozymes = Dict(
    x => Dict(
        "isozyme" => X.IsozymeT{Ex}(
            iso.gene_product_stoichiometry,
            rid_kcat[Symbol(x)],
            nothing,
        ) for (_, iso) in y
    ) for (x, y) in pruned_reaction_isozymes
)

pkm = X.enzyme_constrained_flux_balance_constraints( # kinetic model
    pruned_model;
    reaction_isozymes = parameter_isozymes,
    gene_product_molar_masses,
    capacity,
)

pruned_solution = D.optimized_values(
    pkm,
    parameter_values;
    objective = pkm.objective.value,
    optimizer = T.Optimizer,
)

# All reactions with zero flux have been removed in the pruned model

sort(collect(pruned_solution.tree.fluxes), by = ComposedFunction(abs, last))

# Genes with zero concentration have also been removed.

sort(abs.(collect(values(pruned_solution.tree.gene_product_amounts))))

# ## Optimal flux modes (OFMs)

# We now wish to find the optimal flux modes (OFMs) of this optimal solution.

# In our model, we have a fixed ATP maintenance flux

pruned_model.reactions["ATPM"].lower_bound

# As a result of this fixed flux, standard EFM theory does not hold. Instead,
# we investigate the optimal flux modes (OFMs). These are the flux modes 
# that are able to produce the optimal rate of the objective, whilst still 
# being constrained to the fixed ATPM.

# To find the OFMs, we must augment the stoichiometric matrix, and can then 
# use the same `find_efms` function, since this augmented matrix behaves 
# as required for the algorithm to work.

N = A.stoichiometry(pruned_model)

# Take the fixed fluxes (here only ATPM) and the objective flux

atpm_idx = findfirst(x -> x == "ATPM", A.reactions(pruned_model))
biomass_idx = findfirst(x -> x == "BIOMASS_Ecoli_core_w_GAM", A.reactions(pruned_model))
fixed_fluxes = [atpm_idx, biomass_idx]
flux_values = [
    pruned_solution.tree.fluxes["ATPM"],
    pruned_solution.tree.fluxes["BIOMASS_Ecoli_core_w_GAM"],
]

# Calculate the OFMs, using the `get_ofms` function

OFMs = get_ofms(Matrix(N), fixed_fluxes, flux_values)
OFM_dicts = [
    Dict(A.reactions(pruned_model) .=> OFMs[:, 1]),
    Dict(A.reactions(pruned_model) .=> OFMs[:, 2]),
]

# Scale the OFMs to have one unit flux through biomass

OFM_dicts = [
    Dict(x => y / OFM_dicts[1]["BIOMASS_Ecoli_core_w_GAM"] for (x, y) in OFM_dicts[1]),
    Dict(x => y / OFM_dicts[2]["BIOMASS_Ecoli_core_w_GAM"] for (x, y) in OFM_dicts[2]),
]

# We see that the first OFM is releasing ethanol but not acetate

OFM_dicts[1]["EX_etoh_e"]
OFM_dicts[1]["EX_ac_e"]

@test OFM_dicts[1]["EX_ac_e"] ≈ 0 #src 
@test OFM_dicts[1]["EX_etoh_e"] > 1e-3 #src


# And the second OFM is releasing acetate but no ethanol

OFM_dicts[2]["EX_etoh_e"]
OFM_dicts[2]["EX_ac_e"]

@test OFM_dicts[2]["EX_etoh_e"] ≈ 0 #src
@test OFM_dicts[2]["EX_ac_e"] > 1e-3 #src


# ## Differentiate OFM usage

# The weighted sum of these two OFMs is equal to the whole optimal flux solution,
# and can be easily calculated using two non-shared reactions


M = [
    OFM_dicts[1]["EX_etoh_e"] OFM_dicts[2]["EX_etoh_e"]
    OFM_dicts[1]["EX_ac_e"] OFM_dicts[2]["EX_ac_e"]
]

v = [
    pruned_solution.tree.fluxes["EX_etoh_e"]
    pruned_solution.tree.fluxes["EX_ac_e"]
]

λ = M \ v

# `λ` tells us the weighting of the two OFMs, so here we have 0.6 units of flux
# through OFM₁ and 0.1 units of flux through OFM₂ to produce our optimal biomass.

# It is of interest to see how changes in the kinetic parameters affect this
# optimal weighting.

parameters = Ex.(collect(keys(parameter_values)))
rid_gcounts = Dict(
    rid => [v.gene_product_stoichiometry for (k, v) in d][1] for
    (rid, d) in pruned_reaction_isozymes
)
rid_pid =
    Dict(rid => [iso.kcat_forward for (k, iso) in v][1] for (rid, v) in parameter_isozymes)
sens = differentiate_efm(
    OFM_dicts,
    parameters,
    rid_pid,
    parameter_values,
    rid_gcounts,
    capacity,
    gene_product_molar_masses,
    T.Optimizer,
)

# To get control coefficients instead of sensititivities, we scale these derivatives.

scaled_sens = Matrix(undef, size(sens, 1), size(sens, 2))
for (i, col) in enumerate(eachcol(sens))
    scaled_sens[:, i] = collect(values(parameter_values))[i] .* col ./ λ
end

# We plot the control coefficients to get a visual overview of the system.

sens_perm = sortperm(scaled_sens[1, :])
scaled_sens[:, sens_perm]
parameters[sens_perm]

f, a, hm = heatmap(
    scaled_sens[:, sens_perm]';
    colormap = CairoMakie.Reverse(:RdBu),
    axis = (
        yticks = (1:2, ["Ethanol producing", "Acetate producing"]),
        yticklabelrotation = pi / 2,
        xticklabelrotation = pi / 2,
        xlabel = "Parameters",
        ylabel = "Control coefficient, param/λ * ∂λ/∂param",
        xticks = (1:length(parameters), string.(parameters[sens_perm])),
        title = "Control coefficients of OFM weightings in optimal solution",
        ylabelsize = 23,
        xlabelsize = 23,
        titlesize = 25,
    ),
);
Colorbar(f[:, end+1], hm)
f

# We see that most kinetic parameters have little effect on the optimal OFM weightings,
# see the reactions from RPE to PGI. Those parameters that do affect optimal weightings
# always increase the use of one OFM and decrease the use of the other.

# Increasing the turnover number of PGK is predicted to decrease the optimal flux through
# the acetate producing OFM, and increase the optimal flux through the ethanol producing OFM.
