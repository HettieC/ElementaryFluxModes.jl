"""
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


function rational_nullspace(A::Matrix)
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
    return Z, pivotrows
end


function remove_linearly_dep_rows(A)
    rA = rref!(copy(Array(A)))
    idxs = Int[]
    for (i,row) in enumerate(eachrow(rA))
        if !all(abs.(row) .<= 1e-6) # remove rows of all zero
            push!(idxs, i)
        end
    end
    return rA[idxs, :], idxs
end

"""
Test the adjacency of two extreme rays 'r1', 'r2', where 'K' is the list of already 
checked row indices, 'd' is the total number of row indices to be checked
"""
function fast_adjacency(r1::Vector{Float64},r2::Vector{Float64},K::Vector{Int64},d::Int64)
    zeta = [i for i in 1:length(r1) if r1[i] == 0 && r2[i] == 0]
    if length(intersect(zeta,K)) >= d - 2
        return true
    else
        return false
    end
end

"""
Helper function to reorder the rows of the nullspace so that it is in the form 
[I; K]. If rational==true, the input matrix and output are rational.
"""
function reorder_ns(A::Matrix;rational=false)
    perm_vec = Int64[]
    if !rational
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
    else
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
    end
    return A[perm_vec,:], perm_vec
end


# function makeBitmap(M::Matrix)
#     mask = Array{Bool}(undef,size(M,1),size(M,2))
#     for i in 1:size(M,1) 
#         mask[i,:] = [j == 0 ? false : true for j in M[i,:]]
#     end
#     return mask
# end


# function main_alg(R, q, p)
#     qsplit = q 
#     R1 = makeBitmap(R[1:qsplit-m,:])
#     R2 = R[qsplit-m+1:qsplit,:]
#     numr = qsplit - m
#     for p in q-m+1:q 
#         new_numr = numr
#         jneg = findall(x -> x < 0, R2[1,:])
#         jpos = findall(x -> x > 0, R2[1,:])
#         for k in 1:length(jneg)
#             for l in 1:length(jpos)
#                 newr = or(R1[:,k],R1[:,l]) # element-wise OR 
#                 if numberOfNullBits(newr) + 1 < qsplit-m-1 
#                     continue 
#                 end 

#                 #adjacency test 
#                 adj = true
#                 r = 0 
#                 while adj && r <= numr 
#                     r = r+1 
#                     testr = or(newr, R1[:,r])
#                     if r != l && r != k && all(testr == newr)
#                         adj = false 
#                     end
#                 end

#                 if adj # if adjacent then combine 
#                     new_numr = new_numr + 1 
#                     R1[:,new_numr] = newr 
#                     R2[:,new_numr] = R2[1,l]*R2[:,k] - R2[1,k]*R2[:,l]
#                 end
#             end
#         end
#         #deletion of negative rays 
#         R1[:,jneg] = [] 
#         R2[:,jneg] = [] 
#         numr = new_numr - length(jneg)

#         #transfer current (ie first) row of R2 as bitmap into R1, then delete this row in R2 
#         R1[p,:] = makeBitmap(R2[1,:])
#         R2[1,:] = [] 
#     end
#     return R1
# end


# function preprocessing(NS,pivotrows, optimizer)
#     # create empty E to collect EFMs 
#     E = [] 
#     # create IFF and FAF, binary vectors with r rows, all elements=false 
#     IFF = repeat([false],size(NS,1))
#     FAF = repeat([false],size(NS,1))

#     # set IFF with indices pivotRows[1:end-1] to true 
#     IFF[pivotrows[1:end-1]] .= true

#     # find first EFM in this present IFF config. if any 
#     m = size(NS,1) #note that m==r number of reactions
#     for i in pivotrows[end]:m
#         FAF[i] = true 
#         flux, isFeasible = checkByLP(IFF, FAF, NS, optimizer) 
#         if isFeasible #feasible solution found in current IFF and FAF config 
#             println("IFF:  ", IFF)
#             println("FAF:  ", FAF)
#             push!(E,flux)
#             break
#         end # failed feasibility test 
#         FAF[i] = false
#     end
#     return E, IFF, FAF
# end

# function checkByLP(IFF,FAF,NS,optimizer)
#     m,n = size(NS)
#     # set up LP constraint matrices and empty objective function 
#     A = [NS; NS[FAF,:]]
#     #b = [zeros(m); ones(sum(FAF))]
#     b = [zeros(m); zeros(sum(FAF))]
#     Aeq = NS[IFF,:]
#     beq = zeros(sum(IFF))
#     f = zeros(n) 
#     model = Model(optimizer.Optimizer)
#     @variable(model,t[1:length(f)])
#     @objective(model,Min,dot(f,t))
#     @constraint(model, A*t .>= b .+1e-8 )
#     @constraint(model, Aeq*t .== beq)
#     @constraint(model, t .>= 0)
#     optimize!(model)
#     isFeasible = termination_status(model) == OPTIMAL
#     flux = NS*value.(t)
#     #model = nothing
#     return flux, isFeasible
# end

# """ 
# Run FVA.
# """
# function runFVA(IFF, FAF, NS, optimizer)
#     # first run FBA to get maximum of objective fn 
#     m,n = size(NS)
#     # set up LP constraint matrices and empty objective function
#     A = [NS; NS[FAF,:]]
#     #b = [zeros(m); ones(sum(FAF))]
#     b = [zeros(m); zeros(sum(FAF))]
#     Aeq = NS[IFF,:]
#     beq = zeros(sum(IFF))
#     f = zeros(n) 
#     model = Model(optimizer.Optimizer)
#     @variable(model,t[1:length(f)])
#     # @objective(model,Max,dot(f,t))
#     # @constraint(model, A*t .<= b)
#     @objective(model,Min,dot(f,t))
#     @constraint(model, A*t .>= b .+ 1e-8)
#     @constraint(model, Aeq*t .== beq)
#     @constraint(model, t .>= 0)
#     optimize!(model)
#     # add a test for if the problem is feasible
#     #isFeasible = termination_status(model) == OPTIMAL
#     obj = objective_value(model)
#     @constraint(model, dot(f,t) == obj)

#     min_max = []
#     # loop through all reactions and find their max and min values

#     for i in 1:length(f)
#         # find min flux
#         @objective(model, Min, t[i])
#         optimize!(model)
#         minFlux = objective_value(model)

#         #find max flux
#         @objective(model, Max, t[i])
#         optimize!(model)
#         maxFlux = objective_value(model)
#         min_max = vcat(min_max,[minFlux maxFlux])
#     end
#     return min_max
# end

# """ 
# Turn E into a vector of Vector{Bool} with true if a reaction is used in EFM.
# """
# function cleanE(E::Vector,tol::Float64)
#     newE = Vector{Vector{Bool}}()
#     for x in E 
#         temp = Vector{Bool}()
#         for (i,y) in enumerate(x)
#             if abs(y) > tol 
#                 push!(temp,true) 
#             else
#                 push!(temp,false) 
#             end
#         end
#         push!(newE,temp)
#     end
#     return unique!(newE) 
# end

# """
# Check if the search should be stopped.
# """ 
# function checkTermination(IFF,IFFlast,NS)
#     stopSearch = false 
#     # if IFF not empty, don't stop
#     if any(IFF)
#         return stopSearch
#     end
#     m,n = size(NS)
#     #insufficient downstream reactions or independent reactions 
#     if m-IFFlast < n-1 || rank(NS[IFFlast+1:end,:]) < n-1 
#         stopSearch = true 
#     end 
#     return stopSearch
# end

# """ 
# Run the main programme with backtracking and forward tracking.
# """
# function main_programme(E, IFF, FAF, NS, S, optimizer)
#     m,n = size(NS)
#     IFFlast = findlast(IFF)
#     # Main programme
#     # repeat 1 
#     while true
#     # backtracking: find terminal IFF that can have positive flux    
#     # repeat_2 
#         while true
#             IFFlast = findlast(IFF)
#             if isnothing(IFFlast)
#                 return E 
#             end
#             FAF[IFFlast:end] .= false
#             IFF[IFFlast] = false 
#             FAF[IFFlast] = true # swap terminal IFF to FAF
#             # check if sufficient reactions to reach DoF-1 true elements in IFF
#             stopSearch = checkTermination(IFF, IFFlast, NS) 
#             if stopSearch 
#                 return E # stop algorithm if there will never be sufficient IFF
#             end
#             flux, isFeasible = checkByLP(IFF,FAF,NS,optimizer)
#             # if feasible soln found, proceed to forward-tracking
#             # until_2: feasible soln found for current IFF and FAF configuration
#             if isFeasible
#                 break
#             end
#         end
        
#         #forward tracking: find new downstream IFF 
#         for i in IFFlast+1:m 
#             IFF[i] = true 
#             # do rank test on ns[IFF,:] to see if all rxns marked true in IFF are linearly indep
#             if rank(NS[IFF,:]) == sum(IFF)
#                 flux, isFeasible = checkByLP(IFF, FAF, NS, optimizer)
#                 if isFeasible # forecast and checkpoint
#                     if sum(IFF) == n-1 
#                         push!(E, flux) 
#                         break 
#                     end
#                     ## run forecast and checkpoint:
#                     min_max = runFVA(IFF, FAF, NS, optimizer)
#                     # reactions with permanent positive flux
#                     vpos = repeat([false],m)
#                     for j in 1:size(min_max)[1]
#                         if all(x -> x>0, min_max[j,:])
#                             vpos[j] = true 
#                         end
#                     end
#                     v_zeroable_downstream = repeat([true],m) # zeroable downstream rxns can be set to zero
#                     for (j,v) in enumerate(vpos)
#                         if v 
#                             v_zeroable_downstream[j] = false
#                         end
#                     end
#                     v_zeroable_downstream[1:i] .= false # reactions that can be constrained to zero
#                     # allow forward tracking to continue if all check-points are satisifed:
#                     if sum(v_zeroable_downstream) + sum(IFF) >= m-1 #cond 1
#                         temp_NS = [[i for (i,x) in enumerate(v_zeroable_downstream) if x || IFF[i]],:]
#                         if rank(temp_NS) == m-1 #cond 2
#                             if size(S,2) - rank(S[:,vpos]) == 0 #cond 3
#                                 # allow forward tracking to continue
#                                 continue
#                             else
#                                 break
#                             end
#                         end
#                     end
#                     # from cond 3, EFM found even if number of iff did not reach DoF-1 
#                     if size(S,2) - rank(S[:,vpos]) == 1 
#                         push!(E,flux) # EFM found, store EFM 
#                     end

#                 else # failed feasibility test 
#                     IFF[i] = false
#                 end 
#             else #failed rank test 
#                 IFF[i] = false 
#             end
#         end
#         # until_1: main programme stopped in repeat_2 loop
#     end
#     return E 
# end