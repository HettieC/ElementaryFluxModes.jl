using Revise, COBREXA, Gurobi, ElementaryFluxModes, Serialization
using JSON, RowEchelon, LinearAlgebra, Tulip, Clarabel

model = deserialize("data/models/pmodel_3_bounds")
flux_sol = flux_balance_analysis_dict(model,Tulip.Optimizer)
N = Matrix(stoichiometry(model.inner))
N = remove_linearly_dep_rows(N)[1]
for (i,x) in enumerate(N) 
    if abs(x) < norm(N,Inf)*eps(Float64)
        N[i] = 0
    end
end
fixed_fluxes = [11,16]
flux_values = [round(flux_sol["ATPM"],digits=8),floor(flux_sol["BIOMASS_Ecoli_core_w_GAM"],digits=8)]
N1 = N[:,setdiff(1:size(N,2),fixed_fluxes)]
N2 = N[:,fixed_fluxes]
w = N2*flux_values
N1w = hcat(N1,w)
serialize("data/models/3boundsN1w",N1w)

N1w = deserialize("data/models/3boundsN1w")
ns = rational_nullspace(N1w)[1]
nsrref = rref(ns')
K = Matrix(nsrref')
# make small entries zero
for (i,x) in enumerate(K) 
    if abs(x) < norm(K,Inf)*eps(Float64)
        K[i] = 0
    end
end 

R = BinaryDD(N1w,K)

### post-processing
# for each column of R, find non-zero elements and solve N_.non_zero*r = 0 
# calculate a nullspace for the columns of N corresponding to non-zero elements of column 
E = Matrix(undef,size(R,1),size(R,2))
for (i,r) in enumerate(eachcol(R))
    non_zero = findall(x->x!=0,r)
    flux_ns = rational_nullspace(N1w[:,non_zero])[1]
    mode = zeros(size(R,1))
    for (j,x) in zip(non_zero,flux_ns)
        mode[j] = abs(x) < 1e-10 ? 0 : x
    end
    E[:,i] = mode
end

# put fixed rxns in correct positions
EFMs = zeros(0, size(R,2))
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
push!(idxs,(idxs[end][2]+1,size(E,1)-1))
for (k, (i,j)) in enumerate(idxs)
    EFMs = vcat(EFMs, E[i:j,:])
    k == length(idxs) || (EFMs = vcat(EFMs, E[end,:]'))
end



### put correct number into ATPM and biomass using the Schuster Schuster formula
### using N1w, N2, v1, v2





efm_1 = Dict(x => y for (x,y) in zip(reactions(model.inner),EFMs[:,1]))
efm_2 = Dict(x => y for (x,y) in zip(reactions(model.inner),EFMs[:,2]))
efm_3 = Dict(x => y for (x,y) in zip(reactions(model.inner),EFMs[:,3]))

for (i,e) in enumerate([efm_1, efm_2, efm_3])
    open("data/efms/3emf$(i).json","w") do io 
        JSON.print(io,e)
    end
end