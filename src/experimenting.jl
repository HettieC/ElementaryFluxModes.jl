
ns = round.(rational_nullspace(S)[1])
nsrref = rref(ns')
R = nsrref' # Get a nullspace

## rearrange rows of R to make it in form [I K]'
perm_order = Int64[]
for (i,rowI) in enumerate(eachrow(Float64.(Matrix(I(size(R,2))))))
    for (j,rowR) in enumerate(eachrow(R))
        if rowR == rowI
            push!(perm_order,j)
        end
    end
end
append!(perm_order,setdiff(1:size(R,1),perm_order))
R = R[perm_order,:]



function DDStandard2(A::Matrix)
    ns = round.(rational_nullspace(A)[1])
    if size(ns,2) == 0 
        ns = nullspace(A)
        if size(ns,2) == 0 
            return println("No nullspace exists of A")
        end
    end
    nsrref = rref(ns')
    R = nsrref' # Get a nullspace
    ## rearrange rows of R to make it in form [I K]'
    perm_order = Int64[]
    for (i,rowI) in enumerate(eachrow(Float64.(Matrix(I(size(R,2))))))
        for (j,rowR) in enumerate(eachrow(R))
            if rowR == rowI
                push!(perm_order,j)
            end
        end
    end
    append!(perm_order,setdiff(1:size(R,1),perm_order))
    R = R[perm_order,:]
    d,n = size(R)
    ρ = Int64[] 
    while ρ != collect(1:d) 
        for j in 1:d
            d,n = size(R)
            # partition the columns of R
            tau_pos = [i for i in 1:n if R[j,i] > 0] 
            tau_0 = [h for h in 1:n if R[j,h] == 0]
            tau_neg = [k for k in 1:n if R[j,k] < 0]
            tau_adj = [(i,k) for i in tau_pos for k in tau_neg if adjacency_test(R[:,i],R[:,k],R)]
            Rnew = Array{Float64}(undef,d,0)
            for (i,k) in tau_adj
                p = R[:,i]
                q = R[:,k]
                r_ik = p[j]*q - q[j]*p # adjacent therefore extreme ray
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
            push!(ρ,j)
        end
    end
    return R, perm_order
end 

S_convex = fix_fluxes(S,fixed_fluxes,flux_values)
E = DDStandard2(S_convex)
E = E[:,[i for i in 1:size(E,2) if E[end,i] > 0]]
E_full = clean_DD_result(E, fixed_fluxes, flux_values)

# keep only rays with fixed_flux row entry >0


E = DDStandard2(S_convex)[1]