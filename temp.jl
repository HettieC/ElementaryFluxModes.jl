#################################################################
##### Enzyme constrained model with two capacity constraints #### 
N = Matrix(stoichiometry(model))
N = remove_linearly_dep_rows(N)[1]
#N = rationalize.(N,tol=0.00001)
### fixed reaction only ATPM, reaction 14
fixed_fluxes = [14]
flux_values = [8.39]
N1 = N[:,setdiff(1:size(N,2),fixed_fluxes)]
N2 = N[:,fixed_fluxes]
w = N2*flux_values
N1w = hcat(N1,w)
#N1w = rationalize.(N1w,tol=1e-5)
ns = rational_nullspace(N1w)[1]
nsrref = rref(ns')
R = Matrix(nsrref')
#R = rationalize.(nsrref',tol=1e-8)
#R ./= sum(abs.(R), dims=1) 
for (i,x) in enumerate(R) 
    if abs(x) < norm(R,Inf)*eps(Float64)
        R[i] = 0
    end
end
R, row_order = reorder_ns(R)
d,n = size(R)
ρ = [1,2,3] 
#while ρ != collect(1:d) 
for j in maximum(ρ):d
    println(j)
    d,n = size(R)
#    R ./= sum(abs.(R), dims=1)
    tau_pos = [i for i in 1:n if R[j,i] > norm(R,Inf)*eps(Float64)]
    tau_0 = [h for h in 1:n if abs(R[j,h]) <= norm(R,Inf)*eps(Float64)]
    tau_neg = [k for k in 1:n if R[j,k] < -norm(R,Inf)*eps(Float64)]
    tau_adj = [(i,k) for i in tau_pos for k in tau_neg if adjacency_test(R[:,i],R[:,k],R)]
    Rnew = Array{Float64}(undef,d,0)
    !isempty(tau_adj) && println(tau_adj)
    for (i,k) in tau_adj
        println("i: $i, k: $k")
        p = R[:,i] 
        q = R[:,k] 
        #r_ik = rationalize.(p[j]*q,tol=1e-8) - rationalize.(q[j]*p,tol=1e-8) # +- - -+ 
        r_ik = p[j]*q - q[j]*p
        Rnew = hcat(Rnew,r_ik)
    end
    #Rnew ./= sum(abs.(Rnew), dims=1) 
    Rtemp = Array{Float64}(undef,d,0)
    for i in tau_pos
        Rtemp = hcat(Rtemp,R[:,i])
    end
    for h in tau_0 
        Rtemp = hcat(Rtemp,R[:,h])
    end
    Rtemp = hcat(Rtemp,Rnew)
    R = Rtemp
    for (i,x) in enumerate(R) 
        if abs(x) < norm(R,Inf)*eps(Float64)
            R[i] = 0
        end
    end
    push!(ρ,j)
end

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

efm1 = Dict(split(x,"#";limit=2)[1] => y for (x,y) in zip(reactions(model)[row_order],E[:,1]))
efm2 = Dict(split(x,"#";limit=2)[1] => y for (x,y) in zip(reactions(model)[row_order],E[:,2]))
efm3 = Dict(split(x,"#";limit=2)[1] => y for (x,y) in zip(reactions(model)[row_order],E[:,3]))
efm4 = Dict(split(x,"#";limit=2)[1] => y for (x,y) in zip(reactions(model)[row_order],E[:,4]))

for (i,e) in enumerate([efm1, efm2, efm3, efm4])
    open("data/emf_N1w_gecko$(i).json","w") do io 
        JSON.print(io,e)
    end
end

