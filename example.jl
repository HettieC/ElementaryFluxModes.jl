using ElementaryModes
blah([1 2 3;4 5 6])
using LinearAlgebra
S = [
    1 -1 0 -1 0 0 0 0 0 0 0 0 0 0 ;
    0 1 -1 0 -1 1 0 0 0 0 0 0 0 0 ;
    0 0 0 1 1 -1 -1 -1 0 0 0 0 0 0 ;
    0 0 0 0 0 0 0 1 -1 1 -1 0 0 0 ;
    0 0 0 0 0 0 1 0 0 0 -1 0 0 1 ;
    0 0 0 0 0 0 0 0 0 0 1 -1 0 0 ;
    0 0 0 0 0 0 -1 0 0 0 0 1 -1 0 ;
    0 0 0 0 0 0 0 0 0 0 0 1 0 -1 ;
]


for k in 1:length(jneg)
    for l in 1:length(jpos)
        newr = or(R1[:,k],R1[:,1]) # element-wise OR 
        #check minimum number of zeros 
        if numberOfNullBits(newr)+1 < qsplit-m-1 
            continue
        end

        #adjacency test 
        adj = 1 
        r = 0
        while adj <= numr && r <= numr 
            r = r+1
            testr = or(newr, R1[:,r])
            if r!=1 && r!=k && all(testr == newr)
                adj = 0 
            end
        end

        if adj # if adjacent then combine 
            new_numr = new_numr + 1
            R1[:,new_numr] = newr 
            R2[:,new_numr] = R2[1,1]*R2[:,k] - R2[1,k]*R2[:,1]
        end
    end
end