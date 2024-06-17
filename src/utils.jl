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

function make_all_irreversible(S::Matrix{Float64},reversible::Vector{Int64})
    S_irrev = copy(S)
    for j in reversible
        S_irrev = hcat(S_irrev,-S[:,j]) # add new reversed reactions on right hand side
    end
    return S_irrev 
end

""" 
$(TYPEDSIGNATURES)

Return the EFMs in terms of the original stoichiometric matrix, with reversible reactions.
"""
function reversible_EFMs(E::Matrix{Float64},reversible::Vector{Int64})
    for efm in eachcol(E)
        for (i,r) in enumerate(reversible)
            if efm[end-length(reversible)+i] != 0.0
                efm[r] -= efm[end-length(reversible)+i]
            end
        end
    end
    E = E[1:end-length(reversible),:]
    return E[:,[i for i in 1:size(E,2) if !all(x->x==0,E[:,i])]]
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
function fix_fluxes(S::Matrix{Float64},fixed_fluxes::Vector{Int64},flux_values::Vector{Float64})
    if length(fixed_fluxes) != length(flux_values)
        throw("Number of fixed reactions does not match given number of fixed fluxes")
    end
    N1 = S[:,[i for i in 1:size(S,2) if i ∉ fixed_fluxes]]
    N2 = S[:,[i for i in 1:size(S,2) if i ∈ fixed_fluxes]]
    w = N2*flux_values
    return hcat(N1,w)
end

"""
$(TYPEDSIGNATURES)

Input: 
    E: matrix of EFMs that has come from DDStandard(S_convex)
    fixed_fluxes: same vector of indices of fixed fluxes from the input of fix_fluxes 
Output:
    E_cleaned: 
"""
function clean_DD_result(E::Matrix{Float64},fixed_fluxes::Vector{Int64},flux_values::Vector{Float64})
    # last row of E corresponds to the fixed fluxes 
    fixed_row = E[end,:]
    #E = [[1:end-1],:]
    for (r,v) in zip(fixed_fluxes,flux_values)
        E = [E[1:r-1,:] ; v*fixed_row'; E[r:end,:]]
    end
    rescale_value = flux_values[1] / E[fixed_fluxes[1],findfirst(x->x!=0,E[fixed_fluxes[1]])]
    return E[1:end-1,:]*rescale_value
end

"""
$(TYPEDSIGNATURES)

Helper function to calculate a nullspace of the matrix A, with all rational entries.
"""
function rational_nullspace(A::Matrix; tol=norm(A,Inf)*eps(Float64))
    m,n = size(A)
    R, pivotrows = rref_with_pivots(A)
    r = length(pivotrows)
    nopiv = collect(1:n) 
    nopiv = [x for x in nopiv if x ∉ pivotrows]
    Z = zeros(n, n-r)

    if n > r 
        for (j,k) in enumerate(nopiv) 
            Z[k,:] .= I(n-r)[j,:]
        end
        if r > 0 
            Z[pivotrows,:] = -R[1:r, nopiv]
        end
    end

    for (i,x) in enumerate(R) 
        if abs(x) < tol
            R[i] = 0
        end
    end
    return Z, pivotrows
end

"""
$(TYPEDSIGNATURES)

Helper function to remove linearly dependent rows of the matrix A.
"""
function remove_linearly_dep_rows(A::Matrix{Float64};tol=1e-15)
    rA = rref!(copy(Array(A)))
    idxs = Int[]
    for (i,row) in enumerate(eachrow(rA))
        if !all(abs.(row) .<= tol) # remove rows of all zero
            push!(idxs, i)
        end
    end
    return rA[idxs, :], idxs
end


"""
$(TYPEDSIGNATURES)

Helper function to reorder the rows of the nullspace so that it is in the form 
[I; K]. 
"""
function reorder_ns(A::Matrix)
    perm_vec = Int64[]
    for (i,row) in enumerate(eachrow(A))
        if length(perm_vec) == size(A,2)
            append!(perm_vec,[i for i in 1:size(A,1) if i ∉ perm_vec])
            break
        else
            if row == Matrix(1.0I, size(A,2), size(A,2))[length(perm_vec)+1,:]
                push!(perm_vec,i)
            end
        end
    end
    return A[perm_vec,:], perm_vec
end
