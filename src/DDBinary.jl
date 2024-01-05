"""
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
    return 1 * bitmap
end

"""
Return the zero set of a vector
"""
zero_set(vec) = [i for (i, x) in enumerate(vec) if x == 0]


"""
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


""" 
Input: 
N stoichiometric matrix 
K initial nullspace in the form [I;K*] 
Output:
R binary elementary flux modes 
"""
function DDBinary(N, K)
    R_binary = Matrix(1 * I(size(K, 2)))
    already_pos = size(K, 2)
    R_remaining = K[size(R_binary, 1)+1:end, :]
    R = copy(K)
    for k = 1:size(N, 2) # iterate through every reaction and ensure that the EFMs are non-negative
        k <= already_pos && continue
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
            R = vcat(R_binary, R_remaining) #Q: better to keep separate? better to never separate?
        end
    end
    return R
end
