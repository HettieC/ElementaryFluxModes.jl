# # _E. coli_ core model

using ElementaryFluxModes

using AbstractFBCModels
import AbstractFBCModels.CanonicalModel as CM
using JSONFBCModels
using FastDifferentiation
using DifferentiableMetabolism
using COBREXA
using JSON
using Tulip

# It has been proven that enzyme constrained models with K constraints 
# will use a maximum of K EFMs in their optimal solution (de Groot 2019 (**check))

# Here, we take the _E. coli_ core model with two enzyme constraints, find the optimal 
# solution, and then find the EFMs in the optimal solution.


!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

include("../../test/data_static.jl")

model = convert(CM.Model,load_model("e_coli_core.json"))
# unconstrain glucose
model.reactions["EX_glc__D_e"].lower_bound = -1000.0
# constrain PFL to zero
model.reactions["PFL"].upper_bound = 0.0

rid_kcat = Dict(k => FastDifferentiation.Node(Symbol(k)) for (k,_) in ecoli_core_reaction_kcats)
kcats = Symbol.(keys(ecoli_core_reaction_kcats))

parameter_values = Dict{Symbol, Float64}()

reaction_isozymes = Dict{String,Dict{String,ParameterIsozyme}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
float_reaction_isozymes = Dict{String,Dict{String,COBREXA.Isozyme}}() #src
for (rid,rxn) in model.reactions
    grrs = rxn.gene_association_dnf
    isnothing(grrs) && continue # skip if no grr available
    haskey(ecoli_core_reaction_kcats, rid) || continue # skip if no kcat data available
    for (i, grr) in enumerate(grrs)

        kcat = ecoli_core_reaction_kcats[rid] * 3.6 # change unit to k/h
        parameter_values[Symbol(rid)] = kcat

        d = get!(reaction_isozymes, rid, Dict{String,ParameterIsozyme}())
        d["isozyme_$i"] = ParameterIsozyme(
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = rid_kcat[rid],
            kcat_reverse = rid_kcat[rid],
        )
        d2 = get!(float_reaction_isozymes, rid, Dict{String,COBREXA.Isozyme}()) #src
        d2["isozyme_$i"] = COBREXA.Isozyme( #src
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), #src
            kcat_forward = kcat, #src
            kcat_reverse = kcat, #src
        ) #src
    end
end


gene_product_molar_masses = Dict(k => v for (k, v) in ecoli_core_gene_product_masses)

capacity = [
    (
        "membrane",
        [x for (x,y) in ecoli_core_subcellular_location if y == "membrane"],
        15.0,
    ),
    (
        "cytosol",
        [x for (x,y) in ecoli_core_subcellular_location if y == "cytosol"],
        35.0,
    ),
]

km = build_kinetic_model(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity
)

ec_solution, _, _, _ = optimized_constraints_with_parameters(
    km,
    parameter_values;
    objective = km.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

ec_solution

ec_solution_cobrexa = enzyme_constrained_flux_balance_analysis( #src
    model; #src
    reaction_isozymes = float_reaction_isozymes, #src
    gene_product_molar_masses = ecoli_core_gene_product_masses, #src
    capacity, #src
    optimizer = Tulip.Optimizer, #src
) #src

@test isapprox(ec_solution.objective, ec_solution_cobrexa.objective; atol = TEST_TOLERANCE) #src

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
sort(collect(ec_solution.fluxes), by=ComposedFunction(abs, last))

@test any(values(ec_solution.fluxes) .≈ 0) #src

# And also many inactive gene products. 

sort(collect(ec_solution.gene_product_amounts), by=last)

@test any(isapprox.(values(ec_solution.gene_product_amounts), 0, atol=1e-8)) #src

# Let us start by pruning the model.

flux_zero_tol = 1e-6 # these bounds make a real difference!
gene_zero_tol = 1e-6 

pruned_reaction_isozymes =
    prune_reaction_isozymes(reaction_isozymes, ec_solution, flux_zero_tol)

pruned_parameter_values = Dict(x=>y for (x,y) in parameter_values if haskey(pruned_reaction_isozymes,String(x)))

pruned_gene_product_molar_masses =
    prune_gene_product_molar_masses(gene_product_molar_masses, ec_solution, gene_zero_tol)

pruned_model = prune_model(model, ec_solution, flux_zero_tol, gene_zero_tol)

# Not that for EFM calculations to work, all reactions must be in the forward direction,
# so we reverse any reactions running backwards in the solution 
# reverse any backwards reactions, and reverse the associated isozymes.

for (r, rxn) in pruned_model.reactions
    if rxn.lower_bound < 0 && ec_solution.fluxes[r] < 0
        pruned_model.reactions[r].upper_bound = -rxn.lower_bound
        pruned_model.reactions[r].lower_bound = 0
        pruned_model.reactions[r].stoichiometry = Dict(x => -y for (x, y) in rxn.stoichiometry)
        !haskey(pruned_reaction_isozymes,r) && continue 
        original_isozyme = collect(values(pruned_reaction_isozymes[r]))[1]
        pruned_reaction_isozymes[r] = Dict(
            r => ParameterIsozyme(
                original_isozyme.gene_product_stoichiometry,
                original_isozyme.kcat_reverse,
                nothing
            )
        )
    end
end

# Build the pruned kinetic model 

pkm = build_kinetic_model(
    pruned_model;
    reaction_isozymes = pruned_reaction_isozymes,
    gene_product_molar_masses = pruned_gene_product_molar_masses,
    capacity,
)

# Check that this is indeed a unique solution 
pruned_ec_solution, _, _, _ = optimized_constraints_with_parameters(
    pkm,
    parameter_values;
    objective = pkm.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

# TODO check fluxes and gene amounts same 

# ## Optimal flux modes (OFMs)

# In this model, we have a fixed ATP maintenance flux
pruned_ec_solution.fluxes["ATPM"]

# As a result of this fixed flux, EFM theory does not hold. Instead,
# we investigate the optimal flux modes (OFMs). These are the flux modes 
# that are able to produce the optimal rate of the objective, whilst still 
# being constrained to the fixed ATPM.

# To find the OFMs, we must augment the stoichiometric matrix, and can then 
# use the same `find_efms` function, since this augmented matrix behaves 
# as required for the algorithm to work.

N = AbstractFBCModels.stoichiometry(pruned_model)

# Take the fixed fluxes (here only ATPM) and the objective flux for the matrix `N2`

atpm_idx = findfirst(x -> x == "ATPM", AbstractFBCModels.reactions(pruned_model))
biomass_idx = findfirst(x -> x == "BIOMASS_Ecoli_core_w_GAM", AbstractFBCModels.reactions(pruned_model))
fixed_fluxes = [atpm_idx,biomass_idx]
flux_values = [pruned_ec_solution.fluxes["ATPM"],pruned_ec_solution.fluxes["BIOMASS_Ecoli_core_w_GAM"]]

N1 = N[:, setdiff(1:size(N,2), fixed_fluxes)]
N2 = N[:, fixed_fluxes]
w = N2*flux_values
N = Matrix(hcat(N1,w))

# Calculate the OFMs, using the `get_efms` function

OFMs = get_ofms(Matrix(N), fixed_fluxes, flux_values)
OFM_dicts = [
    Dict(AbstractFBCModels.reactions(pruned_model) .=> OFMs[:,1]),
    Dict(AbstractFBCModels.reactions(pruned_model) .=> OFMs[:,2])
]
# We see that only the first OFM is producing acetate, so we can call this the 
# fermentative OFM.

OFM_dicts[1]["EX_ac_e"]
OFM_dicts[2]["EX_ac_e"]

# The second OFM is consuming 5 times as much oxygen as the first,  so we call 
# this second OFM the respiritory OFM. 

OFM_dicts[1]["EX_o2_e"]
OFM_dicts[2]["EX_o2_e"]

# ## Differentiate OFM usage 

# The weighted sum of these two OFMs is equal to the whole optimal flux solution.
# It is of interest to see how changes in the kinetic parameters affects this 
# optimal weighting.

parameters = FastDifferentiation.Node.(collect(keys(pruned_parameter_values)))
param_vals = collect(values(pruned_parameter_values))
rid_gcounts = Dict(rid => [v.gene_product_stoichiometry for (k, v) in d][1] for (rid, d) in pruned_reaction_isozymes)
rid_pid = Dict(rid => [iso.kcat_forward for (k, iso) in v][1] for (rid, v) in pruned_reaction_isozymes)
sens = differentiate_efm(OFM_dicts, parameters, rid_pid, pruned_parameter_values, rid_gcounts, capacity, gene_product_molar_masses, Tulip.Optimizer)

# To get control coefficients instead of sensititivities, we scale these derivatives.
# First, we must calculate the initial optimal weighting, `lambda` of the OFMs.
# To do so, we must choose fluxes that are different amongst the two OFMs, 
# and since we will be solving a system of linear equations it is best to 
# choose one reaction with zero flux in the first OFM, and the other with 
# zero flux in the second.

M = [
    OFM_dicts[1]["EX_ac_e"] OFM_dicts[2]["EX_ac_e"]
    OFM_dicts[1]["ACALD"] OFM_dicts[2]["ACALD"]
]
v = [
    pruned_ec_solution.fluxes["EX_ac_e"]
    pruned_ec_solution.fluxes["ACALD"]
]
lambda = M\v

scaled_sens = Matrix(undef, size(sens, 1), size(sens, 2))
for (i, col) in enumerate(eachcol(sens))
    scaled_sens[:, i] = param_vals[i] .* col ./ lambda
end

# We plot the control coefficients to get a visual overview of the system.

sens_perm = sortperm(scaled_sens[1,:])
scaled_sens[:,sens_perm]
parameters[sens_perm]

heatmap(parameters[sens_perm],["Fermentative","Respiratory"],)

f, a, hm = heatmap(
    scaled_sens[:,sens_perm]';
    axis = (
        yticks = (1:2, ["Respiratory","Fermentative"]),
        yticklabelrotation = pi / 2,
        xticklabelrotation = pi / 2,
        xlabel = "Parameters",
        ylabel = "Control coefficient, p/λ * ∂λ/∂p",
        xticks = (1:length(parameters), string.(parameters[sens_perm])),
        title = "OFM weighting control coefficients",
    ),
);
Colorbar(f[:,end+1],hm)
f 

# We see that all kinetic parameters increase the use of one OFM and decrease the use 
# of the other.
