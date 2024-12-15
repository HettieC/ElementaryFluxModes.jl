# # _E. coli_ core model

using ElementaryFluxModes

using AbstractFBCModels
import AbstractFBCModels.CanonicalModel as CM
using COBREXA
using JSON

# It has been proven that enzyme constrained models with K constraints 
# will use a maximum of K EFMs in their optimal solution (de Groot 2019 (**check))

# Here, we take the _E. coli_ core model with two enzyme constraints, find the optimal 
# solution, and then find the EFMs in the optimal solution.


!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

include("../../test/data_static.jl")

model = convert(CM.Model,load_model("e_coli_core.json"))
# unconstrain glucose
rids = [x["id"] for x in model.reactions]
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

# ## Calculate EFMs of the pruned model 

# TODO write about OFMs

N = AbstractFBCModels.stoichiometry(pruned_model)