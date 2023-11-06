function prune_model(model::StandardModel, reaction_fluxes::Dict{String,Float64}; make_forward = true, atol = 1e-9, verbose = true)
    pruned_model = StandardModel("pruned_model")

    rxns = Vector{Reaction}()
    mets = Vector{Metabolite}()
    gs = Vector{Gene}()
    mids = String[]
    gids = String[]

    for rid in reactions(model)
        abs(reaction_fluxes[rid]) <= atol && continue

        rxn = deepcopy(model.reactions[rid])
        if reaction_fluxes[rid] > 0
            rxn.lb = max(0, rxn.lb)
        elseif make_forward
            rxn.metabolites = Dict(x => -y for (x,y) in rxn.metabolites)
            rxn.ub = -model.reactions[rid].lb
            rxn.lb = max(0, -model.reactions[rid].ub)
        else
            rxn.ub = min(0, rxn.ub)
        end
        push!(rxns, rxn)

        rs = reaction_stoichiometry(model, rid)
        for mid in keys(rs)
            push!(mids, mid)
        end
        grrs = reaction_gene_association(model, rid)
        isnothing(grrs) && continue
        for grr in grrs
            append!(gids, grr)
        end
    end

    for mid in unique(mids)
        push!(mets, model.metabolites[mid])
    end

    for gid in unique(gids)
        push!(gs, model.genes[gid])
    end

    add_reactions!(pruned_model, rxns)
    add_metabolites!(pruned_model, mets)
    add_genes!(pruned_model, gs)

    #: print some info about process
    if verbose
        rrn = n_reactions(model) - n_reactions(pruned_model)
        @info "Removed $rrn reactions."
    end

    return pruned_model
end