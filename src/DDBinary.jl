"""
$(TYPEDSIGNATURES)

Calculate elementary flux modes of a stoichiometric matrix.
Input:
`N`: the stoichiometric matrix of a network with only forward reactions.
Output:
matrix of the EFMs, each column is a different EFM, rows correspond to the reaction indices of the input `N`.
"""
function get_efms(N::Matrix{Float64}; tol = 1e-15)
    N = remove_linearly_dep_rows_qr(N)
    K = rational_nullspace(Matrix(N))[1]
    # Permute the rows of `K` to be in the form `[I;K*]`
    K, order = reorder_ns(K)

    # The reaction order of `N` must match that of `K`
    N = N[:, order]

    # Run the double description algorithm

    R = DDBinary(N, K)

    # for each column r of R, find non-zero elements and solve N_.non_zero*r = 0
    # calculate a nullspace for the columns of N corresponding to non-zero elements of column
    E = Matrix(undef, size(R, 1), size(R, 2))
    for (i, r) in enumerate(eachcol(R))
        non_zero = findall(x -> x != 0, r)
        flux_ns = rational_nullspace(N[:, non_zero]; tol = 1e-14)[1]
        mode = zeros(size(R, 1))
        for (j, x) in zip(non_zero, flux_ns)
            mode[j] = abs(x) < 1e-14 ? 0 : x
        end
        E[:, i] = mode
    end
    # now put back into original order
    E = E[invperm(order), :]
    return E ./ E[end]
end

"""
$(TYPEDSIGNATURES)
Calculate the optimal flux modes of a pruned optimal solution.
Arguments:
- `N`: the stoichiometric matrix of a pruned model
- `fixed_fluxes`: the indexes of the fixed fluxes
- `flux_values`: the optimal values of the fixed fluxes
Output: matrix of the OFMs, each column is a different OFM, rows correspond to the reaction indices of the input `N`.
"""
function get_ofms(
    N::Matrix{Float64},
    fixed_fluxes::Vector{Int},
    flux_values::Vector{Float64},
)
    if length(fixed_fluxes) != length(flux_values)
        error(
            "Length of fixed fluxes ($(length(fixed_fluxes)) does not match length of flux values ($(length(flux_values)))",
        )
    end
    fixed_order = sortperm(fixed_fluxes)
    fixed_fluxes = fixed_fluxes[fixed_order]
    flux_values = flux_values[fixed_order]
    remaining_rxn_idxs = setdiff(1:size(N, 2), fixed_fluxes)
    N = remove_linearly_dep_rows_qr(N)
    N1 = N[:, setdiff(1:size(N, 2), fixed_fluxes)]
    N2 = N[:, fixed_fluxes]
    w = N2 * flux_values
    N = Matrix(hcat(N1, w))

    K = rational_nullspace(Matrix(N))[1]
    # Permute the rows of `K` to be in the form `[I;K*]`
    K, order = reorder_ns(K)

    # The reaction order of `N` must match that of `K`
    N = N[:, order]

    # Run the double description algorithm

    R = DDBinary(N, K)

    # for each column r of R, find non-zero elements and solve N_.non_zero*r = 0
    # calculate a nullspace for the columns of N corresponding to non-zero elements of column
    E = Matrix(undef, size(R, 1), size(R, 2))
    for (i, r) in enumerate(eachcol(R))
        non_zero = findall(x -> x != 0, r)
        flux_ns = ElementaryFluxModes.rational_nullspace(N[:, non_zero]; tol = 1e-14)[1]
        mode = zeros(size(R, 1))
        for (j, x) in zip(non_zero, flux_ns)
            mode[j] = abs(x) < 1e-14 ? 0 : x
        end
        E[:, i] = mode
    end
    # put back into original order
    E = E[invperm(order), :]
    E = E ./ E[end, :]'

    # put the fixed fluxes back in the right order
    new_order = sortperm(vcat(remaining_rxn_idxs, fixed_fluxes...))
    new_E = vcat(E[1:end-1, :], (E[end, :] * flux_values')')[new_order, :]

    return new_E
end

"""
$(TYPEDSIGNATURES)

Implement the Double Description method in binary form.
The input variables are:
-`N`: the stoichiometric matrix with only forward reactions. If any fluxes are fixed then N has the form [N1 w] where N1 is the stoichiometric matrix for non-fixed fluxes, and w is the columns of fixed fluxes multiplied by their flux values.
-`K`: initial nullspace of N, in the form [I;K*]
Output:
-`R`: binary elementary flux modes
"""
function DDBinary(N, K)
    R_binary = Matrix(I(size(K, 2)))
    already_pos = size(K, 2)
    R_remaining = K[size(R_binary, 1)+1:end, :]
    R = copy(K)
    for k = already_pos+1:size(N, 2) # iterate through every reaction and ensure that the EFMs are non-negative
        if all(x -> x >= 0, R[k, :]) ## non-negativity already satisfied
            R_binary = vcat(R_binary, make_bitmap(R[k, :])')
            R_remaining = R_remaining[2:end, :]
            R = vcat(R_binary, R_remaining)
        else ## remove negative columns and create new rays from
            ## adjacent negative and positive rays
            tau_neg = [i for (i, x) in enumerate(R[k, :]) if x < 0]
            tau_pos = [i for (i, x) in enumerate(R[k, :]) if x > 0]
            R = vcat(R_binary, R_remaining)
            if size(R, 2) <= 3 # if only two or three rays they must all be adjacent
                adj = [(i, j) for i in tau_pos for j in tau_neg]
            else
                adj = [(i, j) for i in tau_pos for j in tau_neg if check_adjacency(i, j, R)]
            end
            new_cols = Matrix(undef, size(N, 2), 0) # create new columns with non-negative combination
            # of one positive and one negative ray
            for (i, j) in adj
                p = copy(R[:, i])
                q = copy(R[:, j])
                new_cols = hcat(new_cols, (p[k]) * q - (q[k]) * p)
            end #Q: is list/matrix comprehension faster?
            R = R[:, setdiff(1:size(R, 2), tau_neg)] # remove negative columns
            R = hcat(R, new_cols) # add the new columns
            R_binary = reduce(hcat, (make_bitmap.(eachcol(R[1:k, :])))) #binarise the 1:k rows of R
            R_remaining = R[k+1:end, :]
            #           R_remaining = R[k+1:end, :]./maximum(abs.(R[k+1,:]))
            R = vcat(R_binary, R_remaining) #Q: better to keep separate? better to never separate?
            #println("k: $k \n size(R): $(size(R)) \n \n")
        end
    end
    return R
end

"""
$(TYPEDSIGNATURES)

Helper function to calculate a nullspace of the matrix A, with all rational entries.
"""
function rational_nullspace(A::Matrix; tol = norm(A, Inf) * eps(Float64))
    m, n = size(A)
    R, pivotrows = rref_with_pivots(A)
    r = length(pivotrows)
    nopiv = collect(1:n)
    nopiv = [x for x in nopiv if x ∉ pivotrows]
    Z = zeros(n, n - r)

    if n > r
        for (j, k) in enumerate(nopiv)
            Z[k, :] .= I(n - r)[j, :]
        end
        if r > 0
            Z[pivotrows, :] = -R[1:r, nopiv]
        end
    end

    for (i, x) in enumerate(R)
        if abs(x) < tol
            R[i] = 0
        end
    end
    return Z, pivotrows
end

"""
$(TYPEDSIGNATURES)

Return a bitmap of the input vector, where entries are 1 if the vector entry is
greater than zero, and zero if the entries are zero.
Return an error if any entries of the vector are less than zero.
"""
function make_bitmap(row)
    bitmap = Bool[]
    for x in row
        if x > 0
            push!(bitmap, true)
        elseif x == 0
            push!(bitmap, false)
        else
            throw(DomainError(x, "argument for `make_bitmap` must be nonnegative"))
        end
    end
    return bitmap
end


"""
$(TYPEDSIGNATURES)

Return the zero set of a vector
"""
zero_set(vec) = [i for (i, x) in enumerate(vec) if x == 0]

"""
$(TYPEDSIGNATURES)

Check adjacency of columns i and j in r by making sure that there exists no other
extreme ray in R whose zero set is a superset of the intersection of the zero sets
of ray i and ray j.
"""
function check_adjacency(i, j, R)
    z1 = zero_set(R[:, i])
    z2 = zero_set(R[:, j])
    z = intersect(z1, z2)
    adj = true
    for col in eachcol(R[:, 1:end.∉[[i, j]]])
        if issubset(z, zero_set(col))
            adj = false
            break
        end
    end
    return adj
end
