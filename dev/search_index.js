var documenterSearchIndex = {"docs":
[{"location":"1-toy-model/","page":"Calculating EFMs in a toy model","title":"Calculating EFMs in a toy model","text":"EditURL = \"1-toy-model.jl\"","category":"page"},{"location":"1-toy-model/#Calculating-EFMs-in-a-toy-model","page":"Calculating EFMs in a toy model","title":"Calculating EFMs in a toy model","text":"","category":"section"},{"location":"1-toy-model/","page":"Calculating EFMs in a toy model","title":"Calculating EFMs in a toy model","text":"using ElementaryFluxModes\n\nusing AbstractFBCModels\nusing LinearAlgebra","category":"page"},{"location":"1-toy-model/#Load-a-simple-model","page":"Calculating EFMs in a toy model","title":"Load a simple model","text":"","category":"section"},{"location":"1-toy-model/","page":"Calculating EFMs in a toy model","title":"Calculating EFMs in a toy model","text":"The code used to construct the model is located in test/simple_model.jl, but it is not shown here for brevity.","category":"page"},{"location":"1-toy-model/","page":"Calculating EFMs in a toy model","title":"Calculating EFMs in a toy model","text":"include(\"../../test/simple_model.jl\"); #hide\n\n\nmodel","category":"page"},{"location":"1-toy-model/","page":"Calculating EFMs in a toy model","title":"Calculating EFMs in a toy model","text":"Get the stoichiometric matrix of the model, this is what we use to find EFMs","category":"page"},{"location":"1-toy-model/","page":"Calculating EFMs in a toy model","title":"Calculating EFMs in a toy model","text":"N = AbstractFBCModels.stoichiometry(model)","category":"page"},{"location":"1-toy-model/#Calculate-the-EFMs","page":"Calculating EFMs in a toy model","title":"Calculate the EFMs","text":"","category":"section"},{"location":"1-toy-model/","page":"Calculating EFMs in a toy model","title":"Calculating EFMs in a toy model","text":"Run the double description algorith on N and K","category":"page"},{"location":"1-toy-model/","page":"Calculating EFMs in a toy model","title":"Calculating EFMs in a toy model","text":"E = get_efms(Matrix(N))","category":"page"},{"location":"1-toy-model/","page":"Calculating EFMs in a toy model","title":"Calculating EFMs in a toy model","text":"If preferred, we can transform the vector of EFMs, E, into a dictionary of reaction => fluxes through the efms","category":"page"},{"location":"1-toy-model/","page":"Calculating EFMs in a toy model","title":"Calculating EFMs in a toy model","text":"EFM_dict = Dict(AbstractFBCModels.reactions(model) .=> eachrow(E))","category":"page"},{"location":"1-toy-model/","page":"Calculating EFMs in a toy model","title":"Calculating EFMs in a toy model","text":"","category":"page"},{"location":"1-toy-model/","page":"Calculating EFMs in a toy model","title":"Calculating EFMs in a toy model","text":"This page was generated using Literate.jl.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"EditURL = \"3-ecoli-core.jl\"","category":"page"},{"location":"3-ecoli-core/#Sensitivities-of-OFMs-in-*E.-coli*-core","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"","category":"section"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"using ElementaryFluxModes\n\nimport AbstractFBCModels as A\nimport AbstractFBCModels.CanonicalModel as CM\nusing JSONFBCModels\nimport FastDifferentiation as F\nconst Ex = F.Node\nimport DifferentiableMetabolism as D\nusing COBREXA\nimport COBREXA as X\nusing JSON\nimport Tulip as T\nusing CairoMakie","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"It has been proven that enzyme constrained models with K constraints will use a maximum of K EFMs in their optimal solution (de Groot 2019).","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"Here, we take the E. coli core model with two enzyme constraints, find the optimal solution, and then find the EFMs in the optimal solution.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"!isfile(\"e_coli_core.json\") &&\n    download(\"http://bigg.ucsd.edu/static/models/e_coli_core.json\", \"e_coli_core.json\")\n\ninclude(\"../../test/data_static.jl\")\n\nmodel = convert(CM.Model, X.load_model(\"e_coli_core.json\"))","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"unconstrain glucose","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"model.reactions[\"EX_glc__D_e\"].lower_bound = -1000.0","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"constrain PFL to zero","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"model.reactions[\"PFL\"].upper_bound = 0.0\n\nrid_kcat = Dict(k => Ex(Symbol(k)) for (k, _) in ecoli_core_reaction_kcats)\nkcats = Symbol.(keys(ecoli_core_reaction_kcats))\n\nfloat_reaction_isozymes = Dict{String,Dict{String,X.Isozyme}}()\nfor (rid, rxn) in model.reactions\n    grrs = rxn.gene_association_dnf\n    isnothing(grrs) && continue # skip if no grr available\n    haskey(ecoli_core_reaction_kcats, rid) || continue # skip if no kcat data available\n    for (i, grr) in enumerate(grrs)\n\n        kcat = ecoli_core_reaction_kcats[rid] * 3.6 # change unit to k/h\n\n        d = get!(float_reaction_isozymes, rid, Dict{String,X.Isozyme}())\n        d[\"isozyme_$i\"] = X.Isozyme(\n            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), # assume subunit stoichiometry of 1 for all isozymes\n            kcat_forward = kcat, # assume forward and reverse have the same kcat\n            kcat_reverse = kcat,\n        )\n    end\nend\n\ngene_product_molar_masses = Dict(k => v for (k, v) in ecoli_core_gene_product_masses)\n\ncapacity = [\n    (\n        \"membrane\", # enzyme group names must be of type Symbol\n        [x for (x, y) in ecoli_core_subcellular_location if y == \"membrane\"],\n        15.0,\n    ),\n    (\"cytosol\", [x for (x, y) in ecoli_core_subcellular_location if y == \"cytosol\"], 35.0),\n]\n\nec_solution = X.enzyme_constrained_flux_balance_analysis(\n    model;\n    reaction_isozymes = float_reaction_isozymes,\n    gene_product_molar_masses,\n    capacity,\n    optimizer = T.Optimizer,\n)\n\nec_solution","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"This solution is both producing acetate and consuming oxygen, therefore it looks like overflow metabolism.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"ec_solution.fluxes[\"EX_ac_e\"]\nec_solution.fluxes[\"EX_o2_e\"]","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"Since we had two enzyme pools, the membrane and the cytosol, this optimal solution will be composed of at most two EFMs. In order to find out how small parameter changes will affect this optimal composition, we must first remove inactive reactions, so that the optimal solution is a unique superposition of the two EFMs, and we can differentiate it.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"This solution contains many inactive reactions","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"sort(collect(ec_solution.fluxes), by = ComposedFunction(abs, last))","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"And also many inactive gene products.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"sort(collect(ec_solution.gene_product_amounts), by = last)","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"Let us start by pruning the model.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"flux_zero_tol = 1e-6 # these bounds make a real difference!\ngene_zero_tol = 1e-6\npruned_model, pruned_reaction_isozymes = D.prune_model(\n    model,\n    ec_solution.fluxes,\n    ec_solution.gene_product_amounts,\n    float_reaction_isozymes,\n    ec_solution.isozyme_forward_amounts,\n    ec_solution.isozyme_reverse_amounts,\n    flux_zero_tol,\n    gene_zero_tol,\n)","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"Since we will later be differentiating our solution, we need to make parameter isozymes and a kinetic model","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"parameter_values = Dict(\n    Symbol(x) => iso.kcat_forward for (x, y) in pruned_reaction_isozymes for (_, iso) in y\n)\n\nrid_kcat = Dict(k => Ex(Symbol(k)) for (k, _) in parameter_values)\nparameter_isozymes = Dict(\n    x => Dict(\n        \"isozyme\" => X.IsozymeT{Ex}(\n            iso.gene_product_stoichiometry,\n            rid_kcat[Symbol(x)],\n            nothing,\n        ) for (_, iso) in y\n    ) for (x, y) in pruned_reaction_isozymes\n)\n\npkm = X.enzyme_constrained_flux_balance_constraints( # kinetic model\n    pruned_model;\n    reaction_isozymes = parameter_isozymes,\n    gene_product_molar_masses,\n    capacity,\n)\n\npruned_solution = D.optimized_values(\n    pkm,\n    parameter_values;\n    objective = pkm.objective.value,\n    optimizer = T.Optimizer,\n)","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"All reactions with zero flux have been removed in the pruned model","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"sort(collect(pruned_solution.tree.fluxes), by = ComposedFunction(abs, last))","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"Genes with zero concentration have also been removed.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"sort(abs.(collect(values(pruned_solution.tree.gene_product_amounts))))","category":"page"},{"location":"3-ecoli-core/#Optimal-flux-modes-(OFMs)","page":"Sensitivities of OFMs in E. coli core","title":"Optimal flux modes (OFMs)","text":"","category":"section"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"We now wish to find the optimal flux modes (OFMs) of this optimal solution.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"In our model, we have a fixed ATP maintenance flux","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"pruned_model.reactions[\"ATPM\"].lower_bound","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"As a result of this fixed flux, standard EFM theory does not hold. Instead, we investigate the optimal flux modes (OFMs). These are the flux modes that are able to produce the optimal rate of the objective, whilst still being constrained to the fixed ATPM.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"To find the OFMs, we must augment the stoichiometric matrix, and can then use the same find_efms function, since this augmented matrix behaves as required for the algorithm to work.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"N = A.stoichiometry(pruned_model)","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"Take the fixed fluxes (here only ATPM) and the objective flux","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"atpm_idx = findfirst(x -> x == \"ATPM\", A.reactions(pruned_model))\nbiomass_idx = findfirst(x -> x == \"BIOMASS_Ecoli_core_w_GAM\", A.reactions(pruned_model))\nfixed_fluxes = [atpm_idx, biomass_idx]\nflux_values = [\n    pruned_solution.tree.fluxes[\"ATPM\"],\n    pruned_solution.tree.fluxes[\"BIOMASS_Ecoli_core_w_GAM\"],\n]","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"Calculate the OFMs, using the get_ofms function","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"OFMs = get_ofms(Matrix(N), fixed_fluxes, flux_values)\nOFM_dicts = [\n    Dict(A.reactions(pruned_model) .=> OFMs[:, 1]),\n    Dict(A.reactions(pruned_model) .=> OFMs[:, 2]),\n]","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"Scale the OFMs to have one unit flux through biomass","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"OFM_dicts = [\n    Dict(x => y / OFM_dicts[1][\"BIOMASS_Ecoli_core_w_GAM\"] for (x, y) in OFM_dicts[1]),\n    Dict(x => y / OFM_dicts[2][\"BIOMASS_Ecoli_core_w_GAM\"] for (x, y) in OFM_dicts[2]),\n]","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"We see that the first OFM is releasing ethanol but not acetate","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"OFM_dicts[1][\"EX_etoh_e\"]\nOFM_dicts[1][\"EX_ac_e\"]","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"And the second OFM is releasing acetate but no ethanol","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"OFM_dicts[2][\"EX_etoh_e\"]\nOFM_dicts[2][\"EX_ac_e\"]","category":"page"},{"location":"3-ecoli-core/#Differentiate-OFM-usage","page":"Sensitivities of OFMs in E. coli core","title":"Differentiate OFM usage","text":"","category":"section"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"The weighted sum of these two OFMs is equal to the whole optimal flux solution, and can be easily calculated using two non-shared reactions","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"M = [\n    OFM_dicts[1][\"EX_etoh_e\"] OFM_dicts[2][\"EX_etoh_e\"]\n    OFM_dicts[1][\"EX_ac_e\"] OFM_dicts[2][\"EX_ac_e\"]\n]\n\nv = [\n    pruned_solution.tree.fluxes[\"EX_etoh_e\"]\n    pruned_solution.tree.fluxes[\"EX_ac_e\"]\n]\n\nλ = M \\ v","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"λ tells us the weighting of the two OFMs, so here we have 0.6 units of flux through OFM₁ and 0.1 units of flux through OFM₂ to produce our optimal biomass.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"It is of interest to see how changes in the kinetic parameters affect this optimal weighting.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"parameters = Ex.(collect(keys(parameter_values)))\nrid_gcounts = Dict(\n    rid => [v.gene_product_stoichiometry for (k, v) in d][1] for\n    (rid, d) in pruned_reaction_isozymes\n)\nrid_pid =\n    Dict(rid => [iso.kcat_forward for (k, iso) in v][1] for (rid, v) in parameter_isozymes)\nsens = differentiate_efm(\n    OFM_dicts,\n    parameters,\n    rid_pid,\n    parameter_values,\n    rid_gcounts,\n    capacity,\n    gene_product_molar_masses,\n    T.Optimizer,\n)","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"To get control coefficients instead of sensititivities, we scale these derivatives.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"scaled_sens = Matrix(undef, size(sens, 1), size(sens, 2))\nfor (i, col) in enumerate(eachcol(sens))\n    scaled_sens[:, i] = collect(values(parameter_values))[i] .* col ./ λ\nend","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"We plot the control coefficients to get a visual overview of the system.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"sens_perm = sortperm(scaled_sens[1, :])\nscaled_sens[:, sens_perm]\nparameters[sens_perm]\n\nf, a, hm = heatmap(\n    scaled_sens[:, sens_perm]';\n    colormap = CairoMakie.Reverse(:RdBu),\n    axis = (\n        yticks = (1:2, [\"Ethanol producing\", \"Acetate producing\"]),\n        yticklabelrotation = pi / 2,\n        xticklabelrotation = pi / 2,\n        xlabel = \"Parameters\",\n        ylabel = \"Control coefficient, param/λ * ∂λ/∂param\",\n        xticks = (1:length(parameters), string.(parameters[sens_perm])),\n        title = \"Control coefficients of OFM weightings in optimal solution\",\n        ylabelsize = 23,\n        xlabelsize = 23,\n        titlesize = 25,\n    ),\n);\nColorbar(f[:, end+1], hm)\nf","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"We see that most kinetic parameters have little effect on the optimal OFM weightings, see the reactions from RPE to PGI. Those parameters that do affect optimal weightings always increase the use of one OFM and decrease the use of the other.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"Increasing the turnover number of PGK is predicted to decrease the optimal flux through the acetate producing OFM, and increase the optimal flux through the ethanol producing OFM.","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"","category":"page"},{"location":"3-ecoli-core/","page":"Sensitivities of OFMs in E. coli core","title":"Sensitivities of OFMs in E. coli core","text":"This page was generated using Literate.jl.","category":"page"},{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/#Binary-DD","page":"Reference","title":"Binary DD","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [ElementaryFluxModes]\nPages = [\"src/DDBinary.jl\"]","category":"page"},{"location":"reference/#ElementaryFluxModes.DDBinary-Tuple{Any, Any}","page":"Reference","title":"ElementaryFluxModes.DDBinary","text":"DDBinary(N, K) -> Any\n\n\nImplement the Double Description method in binary form. The input variables are: -N: the stoichiometric matrix with only forward reactions. If any fluxes are fixed then N has the form [N1 w] where N1 is the stoichiometric matrix for non-fixed fluxes, and w is the columns of fixed fluxes multiplied by their flux values. -K: initial nullspace of N, in the form [I;K*] Output: -R: binary elementary flux modes\n\n\n\n\n\n","category":"method"},{"location":"reference/#ElementaryFluxModes.check_adjacency-Tuple{Any, Any, Any}","page":"Reference","title":"ElementaryFluxModes.check_adjacency","text":"check_adjacency(i, j, R) -> Bool\n\n\nCheck adjacency of columns i and j in r by making sure that there exists no other extreme ray in R whose zero set is a superset of the intersection of the zero sets of ray i and ray j.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ElementaryFluxModes.get_efms-Tuple{Matrix{Float64}}","page":"Reference","title":"ElementaryFluxModes.get_efms","text":"get_efms(N::Matrix{Float64}; tol) -> Any\n\n\nCalculate elementary flux modes of a stoichiometric matrix. Input: N: the stoichiometric matrix of a network with only forward reactions. Output: matrix of the EFMs, each column is a different EFM, rows correspond to the reaction indices of the input N.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ElementaryFluxModes.get_ofms-Tuple{Matrix{Float64}, Vector{Int64}, Vector{Float64}}","page":"Reference","title":"ElementaryFluxModes.get_ofms","text":"get_ofms(\n    N::Matrix{Float64},\n    fixed_fluxes::Vector{Int64},\n    flux_values::Vector{Float64}\n) -> Any\n\n\nCalculate the optimal flux modes of a pruned optimal solution. Arguments:\n\nN: the stoichiometric matrix of a pruned model\nfixed_fluxes: the indexes of the fixed fluxes\nflux_values: the optimal values of the fixed fluxes\n\nOutput: matrix of the OFMs, each column is a different OFM, rows correspond to the reaction indices of the input N.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ElementaryFluxModes.make_bitmap-Tuple{Any}","page":"Reference","title":"ElementaryFluxModes.make_bitmap","text":"make_bitmap(row) -> Vector{Bool}\n\n\nReturn a bitmap of the input vector, where entries are 1 if the vector entry is greater than zero, and zero if the entries are zero. Return an error if any entries of the vector are less than zero.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ElementaryFluxModes.rational_nullspace-Tuple{Matrix}","page":"Reference","title":"ElementaryFluxModes.rational_nullspace","text":"rational_nullspace(\n    A::Matrix;\n    tol\n) -> Tuple{Matrix{Float64}, Any}\n\n\nHelper function to calculate a nullspace of the matrix A, with all rational entries.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ElementaryFluxModes.zero_set-Tuple{Any}","page":"Reference","title":"ElementaryFluxModes.zero_set","text":"zero_set(vec) -> Vector\n\n\nReturn the zero set of a vector\n\n\n\n\n\n","category":"method"},{"location":"reference/#Differentiating","page":"Reference","title":"Differentiating","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [ElementaryFluxModes]\nPages = [\"src/differentiate_efm.jl\"]","category":"page"},{"location":"reference/#ElementaryFluxModes._cost_matrix-Tuple{Vector{Dict{String, Float64}}, Any, Any, Vector{Tuple{String, Vector{String}, Float64}}, Dict{String, Float64}}","page":"Reference","title":"ElementaryFluxModes._cost_matrix","text":"_cost_matrix(\n    EFMs::Vector{Dict{String, Float64}},\n    rid_pid,\n    rid_gcounts,\n    capacity::Vector{Tuple{String, Vector{String}, Float64}},\n    gene_product_molar_masses::Dict{String, Float64};\n    evaluate,\n    parameter_values\n) -> Matrix{Any}\n\n\nCalculate a matrix of the cost vectors of each EFM to each constraint. Entry (i,j) gives the total cost in constraint i to produce one unit objective flux through EFM j. Cost is calculated as ∑w(i)V(j)/kcat, where the variables are:\n\n'w(i)': fraction of the ith enzyme pool that one mole of the enzyme uses up\n'V(j)': flux through the reaction in EFM j\n'kcat': turnover number of the enzyme.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ElementaryFluxModes.differentiate_efm-Tuple{Vector{Dict{String, Float64}}, Vector{FastDifferentiation.Node}, Dict{String, FastDifferentiation.Node}, Dict{Symbol, Float64}, Dict{String, Dict{String, Float64}}, Vector{Tuple{String, Vector{String}, Float64}}, Dict{String, Float64}, Any}","page":"Reference","title":"ElementaryFluxModes.differentiate_efm","text":"differentiate_efm(\n    EFMs::Vector{Dict{String, Float64}},\n    parameters::Vector{FastDifferentiation.Node},\n    rid_pid::Dict{String, FastDifferentiation.Node},\n    parameter_values::Dict{Symbol, Float64},\n    rid_gcounts::Dict{String, Dict{String, Float64}},\n    capacity::Vector{Tuple{String, Vector{String}, Float64}},\n    gene_product_molar_masses::Dict{String, Float64},\n    optimizer\n) -> Any\n\n\nDifferentiate the usage of EFMs with respect to changes in the model parameters using implicit differentiation of the Lagrangian. Calculating the proportion of each EFM used in the optimal solution can be turned into the following LP:     maximise     ∑xᵢ     subject to   Dx = 1 where D is the cost matrix associated to each EFM in each enzyme pool. We then differentiate Lagrangian of this LP to calculate the differential of x by parameters.\n\nArguments:\n\n'EFMs': a vector of dictionaries of EFMs\n'parameters': vector of the model parameters\nrid_pid: dict of reaction ids => parameter ids\nparameter_values: dict of parameter symbol => float value\nrid_gcounts: dict of reaction id => gene product counts\ncapacity: enzyme capacity limit of the enzyme constrained model\ngene_product_molar_masses: dict of gene product geneproductmolar_masses\n\nOutput: Matrix Aᵢⱼ of the the differential of EFM i with respect to parameter j: d(EFMᵢ)/d(pⱼ)\n\n\n\n\n\n","category":"method"},{"location":"reference/#Utilities","page":"Reference","title":"Utilities","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [ElementaryFluxModes]\nPages = [\"src/utils.jl\"]","category":"page"},{"location":"reference/#ElementaryFluxModes.fix_fluxes-Tuple{Matrix{Float64}, Vector{Int64}, Vector{Float64}}","page":"Reference","title":"ElementaryFluxModes.fix_fluxes","text":"fix_fluxes(\n    S::Matrix{Float64},\n    fixed_fluxes::Vector{Int64},\n    flux_values::Vector{Float64}\n) -> Matrix{Float64}\n\n\nFunction to make a convex polyhedron out of a problem with fixed fluxes. Input:     S: irreversible stoichiometric matrix, aka any reversible reactions have         been split into two irreversible reactions     fixedfluxes: list of indices of the fixed fluxes     fluxvalues: list of fixed flux values in same order as fixedfluxes Output:     Sconvex: stoichiometric matrix with which to implement the DD algorithm Note: results of DDStandard need to be transformed to take into account these     fixed fluxes, using cleanDDresult\n\n\n\n\n\n","category":"method"},{"location":"reference/#ElementaryFluxModes.remove_linearly_dep_rows_qr-Tuple{Any}","page":"Reference","title":"ElementaryFluxModes.remove_linearly_dep_rows_qr","text":"remove_linearly_dep_rows_qr(A) -> Any\n\n\nHelper function to remove linearly dependent rows of the matrix A, taken from DifferentiableMetabolism.jl, using QR decomposition.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ElementaryFluxModes.reorder_ns-Tuple{Matrix}","page":"Reference","title":"ElementaryFluxModes.reorder_ns","text":"reorder_ns(A::Matrix) -> Tuple{Matrix, Vector{Int64}}\n\n\nHelper function to reorder the rows of the nullspace so that it is in the form [I; K].\n\n\n\n\n\n","category":"method"},{"location":"reference/#ElementaryFluxModes.reversible_EFMs-Tuple{Matrix{Float64}, Vector{Int64}}","page":"Reference","title":"ElementaryFluxModes.reversible_EFMs","text":"reversible_EFMs(\n    E::Matrix{Float64},\n    reversible::Vector{Int64}\n) -> Matrix{Float64}\n\n\nReturn the EFMs in terms of the original stoichiometric matrix, with reversible reactions.\n\n\n\n\n\n","category":"method"},{"location":"#ElementaryFluxModes.jl","page":"README","title":"ElementaryFluxModes.jl","text":"","category":"section"},{"location":"","page":"README","title":"README","text":"Modules = [ElementaryFluxModes]\nPages = [\"src/ElementaryFluxModes.jl\"]","category":"page"},{"location":"#ElementaryFluxModes.ElementaryFluxModes","page":"README","title":"ElementaryFluxModes.ElementaryFluxModes","text":"Package ElementaryFluxModes provides a Julia implementation of the Double Description method to calculate extreme rays of convex polyhedral cones. We follow the method described in Terzer 2009 thesis for the polyhedral cone Ρ = {x ∈ ℜ^d | Ax = 0, x >= 0}.\n\nThe package can calculate elementary flux modes (EFMs) of optimal solutions to homogeneous enzyme-constrained genome scale metabolic models (ecGSMMs), and optimal flux modes (OFMs) of optimal solutions of inhomogeneous ecGSMMs.\n\nIt is also possible to use differentiation to calculate the sensitivity of the optimal composition of these EFMs or OFMs with respect to model parameters.\n\n\n\n\n\n","category":"module"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"EditURL = \"2-differentiate.jl\"","category":"page"},{"location":"2-differentiate/#Sensitivities-of-the-EFMs-of-a-toy-model","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"","category":"section"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"The optimal flux distribution of any metabolic model can be written as a weighted sum of the EFMs of that model. We are interested in calculating the sensitivity of these weightings to the model parameters, and can use DifferentiableMetabolism.jl to efficiently do so.","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"using ElementaryFluxModes\n\nimport AbstractFBCModels as A\nimport FastDifferentiation as F\nconst Ex = F.Node\nimport DifferentiableMetabolism as D\nusing COBREXA\nimport COBREXA as X\nusing JSON\nimport Tulip as T","category":"page"},{"location":"2-differentiate/#Build-a-simple-enzyme-constrained-model","page":"Sensitivities of the EFMs of a toy model","title":"Build a simple enzyme constrained model","text":"","category":"section"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"The code used to construct the model is located in test/simple_model.jl, but it is not shown here for brevity.","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"include(\"../../test/simple_model.jl\"); #hide\n\nmodel\n\nparameter_values = Dict{Symbol,Float64}()\nreaction_isozymes = Dict{String,Dict{String,X.IsozymeT{Ex}}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.\nfor rid in A.reactions(model)\n    grrs = A.reaction_gene_association_dnf(model, rid)\n    isnothing(grrs) && continue # skip if no grr available\n    haskey(rid_kcat, rid) || continue # skip if no kcat data available\n    for (i, grr) in enumerate(grrs)\n\n        kcat = rid_kcat[rid]\n        parameter_values[Symbol(rid)] = kcat\n\n        d = get!(reaction_isozymes, rid, Dict{String,X.IsozymeT{Ex}}())\n        d[\"isozyme_$i\"] = X.IsozymeT{Ex}(\n            gene_product_stoichiometry = Dict(grr .=> fill(float(1.0), size(grr))), # assume subunit stoichiometry of 1 for all isozymes\n            kcat_forward = Ex(Symbol(rid)),\n            kcat_reverse = 0.0,\n        )\n    end\nend\n\nkm = enzyme_constrained_flux_balance_constraints( # kinetic model\n    model;\n    reaction_isozymes,\n    gene_product_molar_masses,\n    capacity,\n)\n\nec_solution = D.optimized_values(\n    km,\n    parameter_values;\n    objective = km.objective.value,\n    optimizer = T.Optimizer,\n    settings = [X.set_optimizer_attribute(\"IPM_IterationsLimit\", 10_000)],\n)\n\nec_solution.tree.fluxes","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"We have a solution that uses every reaction, and the enzyme capacities are both full. Therefore, we may calculate the EFMs of this solution and directly differentiate them, with no pruning required.","category":"page"},{"location":"2-differentiate/#Calculate-EFMs-of-the-optimal-solution","page":"Sensitivities of the EFMs of a toy model","title":"Calculate EFMs of the optimal solution","text":"","category":"section"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"We need to input the stoichiometric matrix N into ElementaryFluxModes.jl","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"N = A.stoichiometry(model)","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"Calculate a flux matrix of the EFMs, the size of which is (n,k), for n reactions and k EFMs","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"E = get_efms(Matrix(N))","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"Make a dictionary out of the EFM result","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"EFM_dict = Dict(A.reactions(model) .=> eachrow(E))\nEFMs = [\n    Dict(k => v[1] / EFM_dict[\"r6\"][1] for (k, v) in EFM_dict),\n    Dict(k => v[2] / EFM_dict[\"r6\"][2] for (k, v) in EFM_dict),\n]","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"The optimal solution, v, can be written as λ₁EFM₁+λ₂EFM₂=v so that the λ give us the weightings of the two EFMs.","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"Let's calculate λ₁ and λ₂, using reactions r1 and r5, as these are not shared by the EFMs","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"M = [\n    EFM_dict[\"r1\"][1] EFM_dict[\"r1\"][2]\n    EFM_dict[\"r5\"][1] EFM_dict[\"r5\"][2]\n]\n\nv = [\n    ec_solution.tree.fluxes[\"r1\"]\n    ec_solution.tree.fluxes[\"r5\"]\n]\n\nλ = M \\ v","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"The optimal solution is therefore made up of 50 units of flux through EFM₁ and 25 units of flux through EFM₂.","category":"page"},{"location":"2-differentiate/#Differentiate-the-EFMs","page":"Sensitivities of the EFMs of a toy model","title":"Differentiate the EFMs","text":"","category":"section"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"We have calculated the EFMs, and now wish to differentiate their weightings, λ, with respect to the model parameters","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"parameters = Ex.(collect(keys(parameter_values)))\np_vals = collect(values(parameter_values))\nrid_pid =\n    Dict(rid => [iso.kcat_forward for (k, iso) in v][1] for (rid, v) in reaction_isozymes)\n\nsens_efm = differentiate_efm(\n    EFMs,\n    parameters,\n    rid_pid,\n    parameter_values,\n    rid_gcounts,\n    capacity,\n    gene_product_molar_masses,\n    T.Optimizer,\n)","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"","category":"page"},{"location":"2-differentiate/","page":"Sensitivities of the EFMs of a toy model","title":"Sensitivities of the EFMs of a toy model","text":"This page was generated using Literate.jl.","category":"page"}]
}
