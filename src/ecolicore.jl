##### using enzyme constrained ecoli core after pruning
using COBREXA, Gurobi, ElementaryFluxModes, Serialization, RowEchelon
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
    if (bounds(model)[1][i] ∉ [-1000,0,1000]) || (bounds(model)[2][i] ∉ [-1000,0,1000])
        push!(fixed_fluxes,i)
        push!(flux_values,bounds(model)[1][i])
    end
end

## make reduced stoichiometric matrix for reactions in fba solution 
# negative sign if reaction has negative flux
kept_reactions = [(-1)^(fba_sol[r]<1e-6)* i for (i,r) in enumerate(reactions(model)) if abs(fba_sol[r])>1e-6]
rxn_names = [x for x in reactions(model) if abs(fba_sol[x])>1e-6]
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

new_N = remove_linearly_dep_rows(new_N)[1]

## add fixed ATPM 
fixed_new_indices = findall(x->x ∈ fixed_fluxes,abs.(kept_reactions))
N1 = new_N[:,setdiff(1:size(new_N,2),fixed_new_indices)]
N2 = new_N[:,fixed_new_indices]
w = N2 * abs.(flux_values)
N1w = hcat(N1,w)


R = round.(qr(new_N).R; digits = 3)

E = DDStandard(N1w)
N1w_red*E[:,1] == repeat([0],size(N1w_red,1))



### using prune_model 
### using prune_model 
function prune_model(model::StandardModel, reaction_fluxes::Dict{String,Float64}; atol = 1e-9, verbose = true)
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
N1 = N[:,setdiff(1:size(N,2),[10,28])]
N2 = N[:,[10,28]]
v2 = [8.39,10]
N = hcat(N1, N2*v2)

E = DDStandard(N)



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



##### Enzyme constrained model with two capacity constraints 
N = deserialize("data/pmodel_S")
N = remove_linearly_dep_rows(N)[1]

### fixed reaction only ATPM, reaction 14
fixed_fluxes = [14]
flux_values = [8.39]
N1 = N[:,setdiff(1:size(N,2),fixed_fluxes)]
N2 = N[:,fixed_fluxes]
w = N2*flux_values
N1w = hcat(N1,w)

ns = round.(rational_nullspace(N1w)[1],digits = 10)
nsrref = rref(ns')
R = round.(nsrref',digits=10) # Get a nullspace
#R ./= sum(abs.(R), dims=1) 
R = reorder_ns(R)[1]
d,n = size(R)
ρ = [1,2,3] 
#while ρ != collect(1:d) 
for j in maximum(ρ):d
    println(j)
    d,n = size(R)
    R ./= sum(abs.(R), dims=1)
    tau_pos = [i for i in 1:n if R[j,i] > 0]
    tau_0 = [h for h in 1:n if R[j,h] == 0]
    tau_neg = [k for k in 1:n if R[j,k] < 0]
    tau_adj = [(i,k) for i in tau_pos for k in tau_neg if adjacency_test(R[:,i],R[:,k],R)]
    Rnew = Array{Float64}(undef,d,0)
    for (i,k) in tau_adj
        p = R[:,i] # +ve or zero except p[j] +ve
        q = R[:,k] # +ve or zero except q[j] -ve
        r_ik = p[j]*q - q[j]*p # +- - -+ 
        # for (a,x) in enumerate(r_ik)
        #     if x < (1e-18)*maximum(r_ik) ## check if this should be zero
        #         r_ik[a] = 0
        #     end
        # end
        Rnew = hcat(Rnew,r_ik)
    end
    Rnew ./= sum(abs.(Rnew), dims=1) 
    Rtemp = Array{Float64}(undef,d,0)
    for i in tau_pos
        Rtemp = hcat(Rtemp,R[:,i])
    end
    for h in tau_0 
        Rtemp = hcat(Rtemp,R[:,h])
    end
    Rtemp = hcat(Rtemp,Rnew)
    R = Rtemp
    push!(ρ,j)
end



# rescale 
R[:,1] = round.(R[:,1]/R[end,1],digits=12)
R[:,2] = round.(R[:,2]/R[end,2],digits=12)

# put fixed fluxes back into correct position
E = zeros(0, size(R,2))
idxs = Tuple[]
k = 1
for i in fixed_fluxes
    if isempty(idxs) 
        push!(idxs,(1,i-1))
        k += 1
    else
        push!(idxs,(idxs[end][2]+1,i-k))
        k +=1
    end
end
push!(idxs,(idxs[end][2]+1,size(R,1)-1))
for (k, (i,j)) in enumerate(idxs)
    E = vcat(E, R[i:j,:])
    k == length(idxs) || (E = vcat(E, R[end,:]'))
end


N1w*R[:,2]


### visualise these EFMs!
using JSON, COBREXA

model = deserialize("data/pgm")

efm_1_binary = Dict(x => y==0 ? 0 : 1 for (x,y) in zip(reactions(model),E[:,1]))
efm_2_binary = Dict(x => y==0 ? 0 : 1 for (x,y) in zip(reactions(model),E[:,2]))

efm_1 = Dict(x => y for (x,y) in zip(reactions(model),E[:,1]))
efm_2 = Dict(x => y for (x,y) in zip(reactions(model),E[:,2]))

open("data/efm_1.json", "w") do io
    JSON.print(io, efm_1)
end
open("data/efm_2.json", "w") do io
    JSON.print(io, efm_2)
end


open("data/efm_1_binary.json", "w") do io
    JSON.print(io, efm_1_binary)
end
open("data/efm_2_binary.json", "w") do io
    JSON.print(io, efm_2_binary)
end