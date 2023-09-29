using ElementaryFluxModes, COBREXA, Gurobi, LinearAlgebra, RowEchelon

model = load_model(StandardModel,"data/e_coli_core.json")
fba_glc = parsimonious_flux_balance_analysis_dict(model,Gurobi.Optimizer) 
fba_etoh = parsimonious_flux_balance_analysis_dict(
    model,Gurobi.Optimizer;modifications = [
        change_constraint("EX_etoh_e";lb = 1.0,ub=1000.0),
    ]
) # fba with forced ethanol production

glc_model = prune_model(model,fba_glc)
etoh_model = prune_model(model,fba_etoh)

merged_model = StandardModel("Merged respiration and fermentation")
add_reactions!(merged_model,collect(values(glc_model.reactions)))
add_metabolites!(merged_model,collect(values(glc_model.metabolites)))
add_reactions!(merged_model,[r for (k,r) in etoh_model.reactions if k ∉ reactions(merged_model)])
add_metabolites!(merged_model,[m for (k,m) in etoh_model.metabolites if k ∉ metabolites(merged_model)])
merged_model.reactions["EX_etoh_e"].lb = 1.0

fba_sol = parsimonious_flux_balance_analysis_dict(merged_model,Gurobi.Optimizer)
N = Matrix(stoichiometry(merged_model))

N = remove_linearly_dep_rows(N)[1]

#### add fixed fluxes 
findall(x->x=="ATPM",reactions(merged_model))
findall(x->x=="EX_glc__D_e",reactions(merged_model))
findall(x->x=="EX_etoh_e",reactions(merged_model))
fixed_fluxes = [10]#,28]
flux_values = [8.39]#,10]
N1 = N[:,setdiff(1:size(N,2),fixed_fluxes)]
N2 = N[:,fixed_fluxes]
w = N2*flux_values
N1w = hcat(N1,w)
E = DDStandard(N1w)


ns = round.(rational_nullspace(N)[1],digits = 6)
nsrref = rref(ns')
R = nsrref' # Get a nullspace
R ./= sum(abs.(R), dims=1) 
d,n = size(R)
ρ = [1,2,3] 
#while ρ != collect(1:d) 
for j in maximum(ρ):d
    println(j)
    d,n = size(R)
    tau_pos = [i for i in 1:n if R[j,i] > 0]
    tau_0 = [h for h in 1:n if R[j,h] == 0]
    tau_neg = [k for k in 1:n if R[j,k] < 0]
    tau_adj = [(i,k) for i in tau_pos for k in tau_neg if adjacency_test(R[:,i],R[:,k],R)]
    Rnew = Array{Float64}(undef,d,0)
    for (i,k) in tau_adj
        p = R[:,i]
        q = R[:,k]
        r_ik = p[j]*q - q[j]*p ## check if this should be zero
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
R[:,1] = round.(R[:,1]/R[end,1],digits=5)
R[:,2] = round.(R[:,2]/R[end,2],digits=5)

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


efm_1 = Dict(x => y for (x,y) in zip(reactions(merged_model),E[:,1]))
efm_2 = Dict(x => y for (x,y) in zip(reactions(merged_model),E[:,2]))

### visualise these EFMs!
using JSON

open("data/efm_1.json", "w") do io
    JSON.print(io, efm_1)
end
open("data/efm_2.json", "w") do io
    JSON.print(io, efm_2)
end



efm_big = Dict(x => y for (x,y) in zip(reactions(merged_model),R[:,3]))

open("data/efm_big.json","w") do io 
    JSON.print(io,efm_big)
end


function pos_neg(x::Float64)
    if x < -1e-6 
        return -1
    elseif x > 1e-6 
        return 1 
    else
        return 0
    end
end


ns = round.(rational_nullspace(N)[1],digits = 6)
nsrref = rref(ns')
R = nsrref' # Get a nullspace
d,n = size(R)
ρ = [1,2,3] 
#while ρ != collect(1:d) 
for j in maximum(ρ):d
    println(j)
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
    R = pos_neg.(Rtemp)
    push!(ρ,j)
end