using ElementaryFluxModes, LinearAlgebra
N = [
    1 0 0 0 0 -1 -1 -1 0 0 0 0;
    0 1 -1 0 0 1 0 0 -1 1 -1 0;
    0 0 0 0 0 0 1 0 1 -1 0 -1;
    0 0 0 0 0 0 0 1 0 0 0 -1;
    0 0 0 0 -1 0 0 0 0 0 0 1;
    0 0 0 -1 0 0 0 0 0 0 1 1
]
K = rational_nullspace(N)[1]
row_order = [3,6,9,10,11,5,12,8,4,7,1,2]
K = K[row_order,:]

function make_bitmap(row)
    bitmap = Bool[]
    for x in row
        if x > 0
            push!(bitmap,true)
        elseif x == 0 
            push!(bitmap,false)
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
    ### check that no zero set of any other ray in R 
    ### contains z
    adj = true
    for col in eachcol(R[:,1:end .∉ [[i,j]]])
        if issubset(z,zero_set(col))
            adj = false
            break
        end
    end
    return adj
end


## convert upper part of matrix to binary 
R_binary = Matrix(1*I(6))
already_pos = size(R_binary,1)
R_remaining = K[size(R_binary,1)+1:end,:]
R = copy(K)

""" 
Input: 
N stoichiometric matrix 
K initial nullspace in the form [I;K*] 
Output:
R binary elementary flux modes 
"""
function BinaryDDAlgorithm(N,K)
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


### post-processing
# for each column of R, find non-zero elements and solve N_.non_zero*r = 0 
# calculate a nullspace for the columns of N corresponding to non-zero elements of column 
E = Matrix(undef,size(R,1),size(R,2))
for (i,r) in enumerate(eachcol(R))
    non_zero = findall(x->x!=0,r)
    flux_ns = rational_nullspace(N[:,non_zero])[1]
    mode = zeros(size(R,1))
    for (j,x) in zip(non_zero,flux_ns)
        mode[j] = abs(x) < 1e-10 ? 0 : x
    end
    E[:,i] = mode
end

# put fixed rxns in correct positions
EFMs = zeros(0, size(R,2))
idxs = Tuple[]
k = 1
for i in fixed_fluxes
    if isempty(idxs) 
        push!(idxs,(1,i-1))
        k += 1
    else
        push!(idxs,(idxs[end][2]+1,i-k))
        k +=1
    end
end
push!(idxs,(idxs[end][2]+1,size(E,1)-1))
for (k, (i,j)) in enumerate(idxs)
    EFMs = vcat(EFMs, E[i:j,:])
    k == length(idxs) || (EFMs = vcat(EFMs, E[end,:]'))
end

efm_1 = Dict(x => y for (x,y) in zip(reactions(model.inner),EFMs[:,1]))
efm_2 = Dict(x => y for (x,y) in zip(reactions(model.inner),EFMs[:,2]))
efm_3 = Dict(x => y for (x,y) in zip(reactions(model.inner),EFMs[:,3]))

for (i,e) in enumerate([efm_1, efm_2, efm_3])
    open("data/efms/3emf$(i).json","w") do io 
        JSON.print(io,e)
    end
end