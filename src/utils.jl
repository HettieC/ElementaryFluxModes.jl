function initialiseR(S::Matrix)
    m,q = size(S)
    ns = nullspace(S)
    ns_rref, pivotrows = rref_with_pivots(ns')
    R = ns_rref' 
    R = round.(R, digits = 8)
    return R, pivotrows, m, q
end


function main_alg(R, q, p)
    qsplit = q 
    R1 = R[1:qsplit-m,:]
    R2 = R[[qsplit-m:1:qsplit],:]
    numr = qsplit - m
    for p in q-m+1:q 
        new_numr = numr
        jneg = findall(x -> x < 0, R2[1,:])
        jpos = findall(x -> x > 0, R2[1,:])
        for k in 1:length(jneg)
            for l in 1:length(jpos)
                newr = or(R1[:,k],R1[:,l]) # element-wise OR 
                if numberOfNullBits(newr) + 1 < qsplit-m-1 
                    continue 
                end 

                #adjacency test 
                adj = true
                r = 0 
                while adj && r <= numr 
                    r = r+1 
                    testr = or(newr, R1[:,r])
                    if r != l && r != k && all(testr == newr)
                        adj = false 
                    end
                end

                if adj # if adjacent then combine 
                    new_numr = new_numr + 1 
                    R1[:,new_numr] = newr 
                    R2[:,new_numr] = R2[1,l]*R2[:,k] - R2[1,k]*R2[:,l]
                end
            end
        end

    #deletion of negative rays 
    R1[:,jneg] = [] 
    R2[:,jneg] = [] 
    numr = new_numr - length(jneg)

    #transfer current (ie first) row of R2 as bitmap into R1, then delete this row in R2 
    R1[p,:] = makeBitmap(R2[1,:])
    R2[1,:] = [] 
    end
    return R1
end
