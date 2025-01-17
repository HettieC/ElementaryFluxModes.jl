"""
$(TYPEDSIGNATURES)

Function to make all reactions in stoichiometric matrix irreversible, and
optionally prune inactive reactions from a solution and only give the
forward reaction of irreversible reactions in the solution.
Input:
    S: stoichiometric matrix
    reversible: list of reversible reaction indices
Output:
    S_irrev: stoichiometric matrix with reversible reactions split
        into two irreversible reactions
"""

function make_all_irreversible(S::Matrix, reversible::Vector{Int64})
    S_irrev = copy(S)
    for j in reversible
        S_irrev = hcat(S_irrev, -S[:, j]) # add new reversed reactions on right hand side
    end
    return S_irrev
end

"""
$(TYPEDSIGNATURES)

Return the EFMs in terms of the original stoichiometric matrix, with reversible reactions.
"""
function reversible_EFMs(E::Matrix{Float64}, reversible::Vector{Int64})
    for efm in eachcol(E)
        for (i, r) in enumerate(reversible)
            if efm[end-length(reversible)+i] != 0.0
                efm[r] -= efm[end-length(reversible)+i]
            end
        end
    end
    E = E[1:end-length(reversible), :]
    return E[:, [i for i = 1:size(E, 2) if !all(x -> x == 0, E[:, i])]]
end

"""
$(TYPEDSIGNATURES)

Function to make a convex polyhedron out of a problem with fixed fluxes.
Input:
    S: irreversible stoichiometric matrix, aka any reversible reactions have
        been split into two irreversible reactions
    fixed_fluxes: list of indices of the fixed fluxes
    flux_values: list of fixed flux values in same order as fixed_fluxes
Output:
    S_convex: stoichiometric matrix with which to implement the DD algorithm
Note: results of DDStandard need to be transformed to take into account these
    fixed fluxes, using clean_DD_result
"""
function fix_fluxes(
    S::Matrix{Float64},
    fixed_fluxes::Vector{Int64},
    flux_values::Vector{Float64},
)
    if length(fixed_fluxes) != length(flux_values)
        throw("Number of fixed reactions does not match given number of fixed fluxes")
    end
    N1 = S[:, [i for i = 1:size(S, 2) if i ∉ fixed_fluxes]]
    N2 = S[:, [i for i = 1:size(S, 2) if i ∈ fixed_fluxes]]
    w = N2 * flux_values
    return hcat(N1, w)
end


"""
$(TYPEDSIGNATURES)

Helper function to remove linearly dependent rows of the matrix A, taken
from DifferentiableMetabolism.jl, using QR decomposition.
"""
function remove_linearly_dep_rows_qr(A)
    #=
    Filter out linearly dependent constraints using QR decomposition. Since the
    problem solved, assume there are no contradictory constraints.

    See: https://math.stackexchange.com/questions/748500/how-to-find-linearly-independent-columns-in-a-matrix

    Make use of the fact that sparse QR returns a staircase profile with column
    ordering by default. The default tolerance of what is a 0 seems okay to rely
    on. Could be a source of bugs though...
    =#
    Is, Js, Vs = SparseArrays.findnz(sparse(A))

    a = SparseArrays.sparse(Js, Is, Vs)

    t = LinearAlgebra.qr(a)  # do transpose here for QR
    max_lin_indep_columns = 0
    for i in axes(t.R, 2) # depends on preordered QR!
        Is, _ = SparseArrays.findnz(t.R[:, i])
        if isempty(Is) || maximum(Is) != i
            break
        else
            max_lin_indep_columns = i
        end
    end

    return A[sort(t.pcol[1:max_lin_indep_columns]), :] # undo permumation
end

"""
$(TYPEDSIGNATURES)

Helper function to reorder the rows of the nullspace so that it is in the form
[I; K].
"""
function reorder_ns(A::Matrix)
    j = 1
    perm_vec = Int64[]
    for (i, row) in enumerate(eachrow(A))
        j > size(A, 2) && break
        if row == Array(Float64.(I(size(A, 2))))[j, :]
            push!(perm_vec, i)
            j += 1
        end
    end
    append!(perm_vec, [i for i = 1:size(A, 1) if i ∉ perm_vec])

    return A[perm_vec, :], perm_vec
end
