### Implementation of the double description method.
### Input: Matrix A, defining a polyhedral cone Ρ = {x ∈ ℜ^d | Ax = 0, x >= 0}
### Output: Matrix R, extreme rays of Ρ = {x ∈ ℜ^d | x = Rc, c >= 0}


### To do: Bit mask the rays for speed
"""
Function to implement the standard DD algorithm, as described in Terzer 2009 thesis, to find 
extreme rays of a polyhedral cone Ρ = {x ∈ ℜ^d | Ax = 0, x >= 0}. 
"""
function DDStandard(A::Matrix)
    ns = round.(rational_nullspace(A)[1])
    if size(ns,2) == 0 
        ns = nullspace(A)
        if size(ns,2) == 0 
            return println("No nullspace exists of A")
        end
    end
    nsrref = rref(ns')
    R = nsrref' # Get a nullspace
    d,n = size(R)
    ρ = Int64[] 
    while ρ != collect(1:d) 
        for j in 1:d
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
            push!(ρ,j)
        end
    end
    return R
end 
 
function DDStandardFixedFluxes(A::Matrix,fixed_fluxes::Vector{Int64})
    ns = round.(rational_nullspace(A)[1])
    if size(ns,2) == 0 
        ns = nullspace(A)
        if size(ns,2) == 0 
            return println("No nullspace exists of A")
        end
    end
    nsrref = rref(ns')
    R = nsrref' # Get a nullspace
    d,n = size(R)
    ρ = fixed_fluxes
    while ρ != collect(1:d) 
        for j in setdiff(1:d,ρ)
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
            push!(ρ,j)
        end
    end
    return R
end 




""" 
Test if two extreme rays r1 and r2 are adjacent in A, adjacency test from Algorithms in Bioinformatics, p335, Terzer, Stelling, 2006
"""
function adjacency_test(r1::Vector{Float64},r2::Vector{Float64},A::Matrix{Float64})
    zeta = [i for i in 1:length(r1) if r1[i] == 0 && r2[i] == 0]
    if rank(A[zeta,:]) == rank(A)-2
        return true 
    else 
        return false
    end
end


