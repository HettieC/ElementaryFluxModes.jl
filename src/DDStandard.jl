"""
$(TYPEDSIGNATURES)

Function to implement the standard DD algorithm, as described in Terzer 2009 thesis, to find 
extreme rays of a polyhedral cone Ρ = {x ∈ ℜ^d | Ax = 0, x >= 0}. 
Input
    A: the stoichiometric matrix with only forward reactions. If any fluxes are fixed 
        then A has the form [N1 w] where N1 is the stoichiometric matrix for non-fixed fluxes,
        and w is the columns of fixed fluxes multiplied by their flux values.
    Output
    R: an n x K matrix containing K EFMs, fixed fluxes need to be put into the correct 
        position in postprocessing.
    row_order: permutation vector of the rows of the nullspace, matching reactions names to 
        fluxes needs to be done using this row_order
"""
function DDStandard(A::Matrix)
    ns = round.(rational_nullspace(A)[1], digits = 6)
    if size(ns, 2) == 0
        ns = nullspace(A)
        if size(ns, 2) == 0
            return println("No nullspace exists of A; problem has no solution.")
        end
    end
    nsrref = rref(ns')
    R = Matrix(nsrref')
    R, row_order = reorder_ns(R) # put I matrix at start of R
    d, n = size(R)
    for j = n+1:d # first n rows of R are guaranteed >=0 from reordering 
        d, n = size(R)
        tau_pos = [i for i = 1:n if R[j, i] > 0]
        # can be sped up if all cols are pos or zero
        tau_0 = [h for h = 1:n if R[j, h] == 0]
        tau_neg = [k for k = 1:n if R[j, k] < 0]
        tau_adj = [
            (i, k) for i in tau_pos for
            k in tau_neg if rank_adjacency_test(R[:, i], R[:, k], R)
        ]
        Rnew = Array{Float64}(undef, d, 0)
        for (i, k) in tau_adj
            p = R[:, i]
            q = R[:, k]
            r_ik = p[j] * q - q[j] * p
            Rnew = hcat(Rnew, r_ik)
        end
        Rtemp = Array{Float64}(undef, d, 0)
        for i in tau_pos
            Rtemp = hcat(Rtemp, R[:, i])
        end
        for h in tau_0
            Rtemp = hcat(Rtemp, R[:, h])
        end
        Rtemp = hcat(Rtemp, Rnew)
        R = Rtemp
    end
    return R, row_order
end


""" 
$(TYPEDSIGNATURES)

Check if two extreme rays r1 and r2 are adjacent in A, using the rank adjacency test from Algorithms in Bioinformatics, p335, Terzer, Stelling, 2006
"""
function rank_adjacency_test(
    r1::Union{Vector{Float64},Vector{Int64}},
    r2::Union{Vector{Float64},Vector{Int64}},
    A::Union{Matrix{Float64},Matrix{Int64}},
)
    zeta = [i for i = 1:length(r1) if r1[i] == 0 && r2[i] == 0]
    if rank(A[zeta, :]) == rank(A) - 2
        return true
    else
        return false
    end
end
