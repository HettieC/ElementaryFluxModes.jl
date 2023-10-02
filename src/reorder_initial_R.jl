## example: ecoli core pgm model

N = deserialize("data/pmodel_S")
N = remove_linearly_dep_rows(N)[1]

### fixed reaction only ATPM, reaction 14
fixed_fluxes = [14]
flux_values = [8.39]
N1 = N[:,setdiff(1:size(N,2),fixed_fluxes)]
N2 = N[:,fixed_fluxes]
w = N2*flux_values
N1w = hcat(N1,w)

ns = round.(rational_nullspace(N1w)[1],digits = 6)
nsrref = rref(ns')
R = round.(nsrref',digits=10) # Get a nullspace
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
        p = R[:,i]
        q = R[:,k]
        r_ik = p[j]*q - q[j]*p
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