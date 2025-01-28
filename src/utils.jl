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
