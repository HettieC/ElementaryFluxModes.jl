##### using ecoli core 
using COBREXA, Gurobi

model = load_model(StandardModel,"data/e_coli_core.json")
fba_sol = parsimonious_flux_balance_analysis_dict(model,Gurobi.Optimizer)

## get list of reversible reaction indices
reversible = Int64[]
for (i,r) in enumerate(reactions(model))
    if bounds(model)[1][i] < 0 && bounds(model)[2][i] > 0 
        push!(reversible,i)
    end
end

## get list of fixed fluxes 
fixed_fluxes = Int64[]
flux_values = Float64[]
for (i,r) in enumerate(reactions(model))
    if ( bounds(model)[1][i] > 0 && bounds(model)[2][i] > 0 ) || ( bounds(model)[1][i] < 0 && bounds(model)[2][i] < 0 )
        push!(fixed_fluxes,i)
        push!(flux_values,bounds(model)[1][i])
    end
end

## make reduced stoichiometric matrix for reactions in fba solution 
# negative sign if reaction has negative flux
kept_reactions = [(-1)^(fba_sol[r]<1e-6)* i for (i,r) in enumerate(reactions(model)) if abs(fba_sol[r])>1e-6]

## make reduced stoichiometric matrix
N = Matrix(stoichiometry(model))
new_N = N[:,abs.(kept_reactions)]
for (i,j) in enumerate(kept_reactions)
    if j < 0 
        new_N[:,i] = -new_N[:,i]
    end
end

# remove rows for metabolites not involved in any reactions
removed_mets = Int64[]
for (i, row) in enumerate(eachrow(new_N))
    if round.(row) == repeat([0],size(new_N,2))
        push!(removed_mets,i)
    end
end
new_N = new_N[setdiff(1:size(new_N,1),removed_mets),:]

# # remove linearly dependent rows with QR decomp
# R = qr(new_N).R
# idxs = Int[]
# for i = 1:size(R, 1)
#     if !all(abs.(R[i, :]) .<= 1e-6)# remove rows of all zero
#         push!(idxs,i)
#     end
# end
# qr_N = new_N[idxs, :]


## add fixed ATPM 
fixed_new_indices = findall(x->x ∈ fixed_fluxes,kept_reactions)
N1 = new_N[:,setdiff(1:size(new_N,2),fixed_new_indices)]
N2 = new_N[:,fixed_new_indices]
w = N2 * flux_values
N1w = hcat(N1,w)

N_ind, ids = remove_linearly_dep_rows(new_N)

R = round.(qr(new_N).R; digits = 3)

E = DDStandard(N_red)
N1w_red*E[:,1] == repeat([0],size(N1w_red,1))



### using prune_model 
function prune_model(model::StandardModel, reaction_fluxes; atol = 1e-9, verbose = true)
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
        else
            rxn.metabolites = Dict(x => -y for (x,y) in rxn.metabolites)
            rxn.ub = -model.reactions[rid].lb
            rxn.lb = max(0, -model.reactions[rid].ub)
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

pmodel = prune_model(model,fba_sol;atol=1e-6)
flux_balance_analysis_dict(pmodel,Gurobi.Optimizer)

N = Matrix(stoichiometry(pmodel))

## remove linearly dependent rows 
N = remove_linearly_dep_rows(N)[1]
N1 = N[:,setdiff(1:size(N,2),[10])]
N2 = N[:,10]
v2 = [8.39]'
N = hcat(N1, N2*v2)

E = DDStandard(N)


ns = (rational_nullspace(N)[1])
nsrref = rref(ns')
R = Matrix(nsrref') # Get a nullspace
d,n = size(R)
ρ = [1,2] 
while ρ != collect(1:d)
    for j in 3:d
        d,n = size(R)
        tau_pos = [i for i in 1:n if R[j,i] > 0]
        tau_0 = [h for h in 1:n if R[j,h] == 0]
        tau_neg = [k for k in 1:n if R[j,k] < 0]
        tau_adj = [(i,k) for i in tau_pos for k in tau_neg if adjacency_test(R[:,i],R[:,k],R)]
        Rnew = Array{Float64}(undef,d,0)
        for (i,k) in tau_adj
            p = R[:,i]
            q = R[:,k]
            r_ik = p[j]*q - q[j]*p
            Rnew = hcat(Rnew,r_ik)
        end
        Rtemp = Array{Float64}(undef,d,0)
        for i in tau_pos
            Rtemp = hcat(Rtemp,R[:,i])
        end
        for h in tau_0 
            Rtemp = hcat(Rtemp,R[:,h])
        end
        Rtemp = hcat(Rtemp,Rnew)
        R = Rtemp
        push!(ρ,j);
    end
end






#### swap rows in R and columns in N so that R has form [I; K] 
new_order = Int64[]
for (i,row) in enumerate(eachrow(R))
    if length(new_order) == size(R,2)
        append!(new_order,i:size(R,1))
        break 
    end
    if row == Matrix(I(size(R,2)))[i]
        push!(new_order,i)
    end
end 




### remove the futile cycle from model 
remove_reaction!(pmodel,reactions(pmodel)[21])
fba_sol = flux_balance_analysis_dict(pmodel,Tulip.Optimizer)

fba_sol[reactions(pmodel)[35]]