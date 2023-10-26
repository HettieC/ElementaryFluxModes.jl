"""
Return a bitmap of the input vector, where entries are 1 if the vector entry is 
greater than zero, and zero if the entries are zero. 
Return an error if any entries of the vector are less than zero.
"""
function make_bitmap(row)
    bitmap = Bool[]
    for x in row
        if x > 0
            push!(bitmap,true)
        elseif x == 0 
            push!(bitmap,false)
        else 
            throw(DomainError(x, "argument for `make_bitmap` must be nonnegative"))
        end
    end
    return 1*bitmap 
end
zero_set(vec) = [i for (i,x) in enumerate(vec) if x==0] # make bit set zero set


"""
Check adjacency of columns i and j in r
"""
function check_adjacency(i,j,R)
    z1 = zero_set(R[:,i])
    z2 = zero_set(R[:,j])
    z = intersect(z1,z2)
    ### check that no zero set of any other ray in R contains z
    adj = true
    for col in eachcol(R[:,1:end .∉ [[i,j]]])
        if issubset(z,zero_set(col))
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
function DDBinary(N,K)
    R_binary = Matrix(1*I(size(K,2)))
    already_pos = size(K,2)
    R_remaining = K[size(R_binary,1)+1:end,:]
    R = copy(K)
    for k in 1:size(N,2)
        k <= already_pos && continue
        if all(x->x>=0,R[k,:]) ## non-negativity already satisfied
            R_binary = vcat(R_binary,make_bitmap(R[k,:])')
            R_remaining = R_remaining[2:end,:]
            R = vcat(R_binary,R_remaining)
        else ## remove negative columns and create new rays from
            ## adjacent negative and positive rays
            neg = [i for (i,x) in enumerate(R[k,:]) if x<0]
            pos = [i for (i,x) in enumerate(R[k,:]) if x>0]
            R = vcat(R_binary,R_remaining)
            #### if only three vectors then all must be adjacte, speed up!
            adj = Tuple[]
            for i in pos 
                for j in neg
                    if check_adjacency(i,j,R)
                        push!(adj,(i,j))
                    end
                end
            end
            new_cols = Matrix(undef,size(N,2),0)
            for (i,j) in adj 
                new_cols = hcat(
                    new_cols, 
                    (R[:,i][k])*R[:,j] - (R[:,j][k])*R[:,i]
                )
                #new_cols = hcat(new_cols,[x+y for (x,y) in zip(R[:,i],R[:,j])]) ### combining maybe needs to be same as original alg?
            end
            R = R[:,setdiff(1:size(R,2),neg)]
            R = hcat(R,new_cols)
            R_binary = reduce(hcat,(make_bitmap.(eachcol(R[1:k,:]))))
            R_remaining = R[k+1:end,:]
            R = vcat(R_binary,R_remaining)
        end
    end
    return R 
end