### Implementation of the double description method.
### Input: Matrix A, defining a polyhedral cone Ρ = {x ∈ ℜ^d | Ax = 0, x >= 0}
### Output: Matrix R, extreme rays of Ρ = {x ∈ ℜ^d | x = Rc, c >= 0}


### To do: Bit mask the rays for speed
"""
Function to implement the standard DD algorithm, as described in Terzer 2009 thesis, to find 
extreme rays of a polyhedral cone Ρ = {x ∈ ℜ^d | Ax = 0, x >= 0}. 
Input
    A: the stoichiometric matrix with only forward reactions, so if any fluxes are fixed 
        then A has the form [N1 w] where N1 is the stoichiometric matrix for non-fixed fluxes,
        and w is the columns of fixed fluxes multiplied by their flux values.#
    row_order: permutation vector of the rows of the nullspace, matching reactions names to 
        fluxes needs to be done using this row_order
Output
    R: an n x K matrix containing K EFMs, fixed fluxes need to be put into the correct 
        position in postprocessing.
"""
function DDStandard(A::Matrix)
    ns = round.(rational_nullspace(A)[1],digits = 6)
    if size(ns,2) == 0 
        ns = nullspace(A)
        if size(ns,2) == 0 
            return println("No nullspace exists of A; problem has no solution.")
        end
    end
    nsrref = rref(ns')
    R = Matrix(nsrref') 
    R, row_order = reorder_ns(R) # put I matrix at start of R
    d,n = size(R)
    for j in n+1:d # first n rows of R are guaranteed >=0 from reordering 
        # println("$(j)th iteration. \n row j-1")
        # println(round.(R[j,:],digits=5))
        d,n = size(R)
        tau_pos = [i for i in 1:n if R[j,i] > 0]
        # can be sped up if all cols are pos or zero
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
        # println("row $j after iteration")
        # println(round.(R[j,:],digits=5))
    end
    return R, row_order
end 

"""
Testing iteration by iteration
"""
function DDStandardIteration(A::Matrix,J::Int64)
    ns = round.(rational_nullspace(A)[1],digits = 6)
    if size(ns,2) == 0 
        ns = nullspace(A)
        if size(ns,2) == 0 
            return println("No nullspace exists of A; problem has no solution.")
        end
    end
    nsrref = rref(ns')
    R = Matrix(nsrref') 
    R, row_order = reorder_ns(R) # put I matrix at start of R
    d,n = size(R)
    for j in n+1:J # first n rows of R are guaranteed >=0 from reordering 
        # println("$(j)th iteration. \n row j-1")
        # println(round.(R[j,:],digits=5))
        d,n = size(R)
        tau_pos = [i for i in 1:n if R[j,i] > 0]
        # can be sped up if all cols are pos or zero
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
        # println("row $j after iteration")
        # println(round.(R[j,:],digits=5))
    end
    return R, row_order
end 

""" 
Test if two extreme rays r1 and r2 are adjacent in A, adjacency test from Algorithms in Bioinformatics, p335, Terzer, Stelling, 2006
"""
function adjacency_test(r1::Union{Vector{Float64},Vector{Int64}},r2::Union{Vector{Float64},Vector{Int64}},A::Union{Matrix{Float64},Matrix{Int64}})
    zeta = [i for i in 1:length(r1) if r1[i] == 0 && r2[i] == 0]
    if rank(A[zeta,:]) == rank(A)-2
        return true 
    else 
        return false
    end
end

"""
Test if two rays r1 and r2 are adjacent in A, test from Terzer 2009, Propn 2.6
"""
function adjacency_test_2(r1,r2,A)
    zeta = [i for (i,x) in enumerate(r1) if x == 0 && r2[i] == 0]
    if rank(hcat(A,Matrix(1.0I(size(R,2)))[zeta,:])) == size(A,1) - 2
        return true
    else
        return false
    end
end