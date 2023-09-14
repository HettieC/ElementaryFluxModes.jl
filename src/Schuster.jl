### Algorithm for when some fluxes are fixed ###

"""
Input: 
    N: full stoichiometric matrix
    v2: fixed fluxes
    N2: stoichiometry of reactions with fixed fluxes 
    N1w: N1*(N2*v2) where N1 is the stoichiometric matrix with fixed fluxes removed 
"""

# Initially, do algorithm with no fixed fluxes, aka w = 0
stoich = copy(S)
nmetab, nflux = size(stoich)
ntotal = nmetab + nflux 
ntab = nflux*2
m = nflux + 1
# say that r1 is a fixed flux 
N2 = stoich[:,1]
v2 = [1]
w = N2*v2'


tol = 1e-10 #limit for non-zero


# intitialisation of Tab
Tab = Matrix(undef,nflux,ntotal)
for i in 1:nflux
    for j in 1:ntotal
        if j > nflux 
            Tab[i,j] = stoich[j-nflux,i]
        elseif i != j 
            Tab[i,j] = 0 
        else 
            Tab[i,j] = 1
        end
    end
end
n_zero = 0 # number of column elmts being zero 
comb = nflux # number of pairs of col vectors with opposite sign


## Own implementation of algorithm 



























B1 = []
B2 = []
C = []
Z = []


### my own implementation of the algorithm 
for i in 1:nmetab
    j = 0
    ntabold = n_zero + comb 
    comb = n_zero 
    breaktest = 0 
    while j < ntabold
        j += 1
        if abs(Tab[j,nflux + i]) < tol # find zero elmts of nflux+ith column of Tab
            for k in 1:ntotal 
                Z[n_zero,k] = Tab[j,k]
                n_zero += 1
            end
        else 
            for k in j+1:ntabold # search for elmts with opposite sign 
                sgn = 0 
                if abs(Tab[j,nflux+i]) > tol && abs(Tab[k,nflux+i]) > tol
                    sgn = Tab[j,nflux+i]*Tab[k,nflux+i]
                end 
                if sgn < 0 
                    help1 = 0 
                    for l in 1:nflux+i 
                        if abs(Tab[j,l]) < tol && abs(Tab[k,l]) < tol 
                            help1 += 1
                            C[help1] = 1
                        end
                    end
                    for l in 1:ntabold # test for condition (14) 
                        help2 = 0
                        if l != j && l != k 
                            for m in 1:help1
                                if abs(Tab[l,C[m]]) > tol 
                                    break 
                                else help2 +=1 
                                end
                            end
                            if help2 == help1 
                                l = ntabold - 1 
                            end
                        end
                        if l == ntabold-1
                            break
                        end 
                        if help1 != help2 
                            B1[comb] = j 
                            comb += 1 
                            B2[comb] = k # cond (14) fulfilled, store (j,k) in B1 and B2 
                        end
                    end
                end
            end 
        end
    end 
    if comb + n_zero > ntab # if number of rows of new matrix Tab is greater than ntab,
        # ntab will be increased and the algorithm starts again 
        B1 = []
        B2 = []
        Tab = []
        Z = []
        ntab += nflux 
        breaktest = 1 
        break 
    end

    for j in 1:comb 
        help1 = B1[j]
        help2 = B2[j] 
        sgn = abs(Tab[help2,nflux+i])>abs(Tab[help1,nflux+i]) ? abs(Tab[help2,nflux+i]) : abs(Tab[help1,nflux+i])
    end
end





















# calculation of subsequent Tab
for i in 1:nmetab
    j = -1 # maybe 0?
    ntabold = comb + n_zero 
    comb = n_zero # 
    breaktest = 0
    while j < ntabold 
        j += 1 
        if abs(Tab[j,nflux+i]) < EPS # find non-zero elmts of (nflux+i)th col of Tab
            for k in 1:ntotal 
                Z[n_zero,k] = Tab[j,k] # store corresponding rows in Z
                n_zero += 1
            end
        else 
            for k in j+1:ntabold # search for elmts with opposite sign 
                sgn = 0 
                if abs(Tab[j,nflux+i]) > EPS && abs(Tab[k,nflux+i]) > EPS 
                    sgn = Tab[j,nflux+i]*Tab[k,nflux+i]
                end
                if sgn < 0 
                    help1 = 0 
                    for l in 1:nflux+i 
                        if abs(Tab[j,l]) < EPS && abs(Tab[k,k]) < EPS
                            C[help1+1] = 1
                        end
                    end
                    for l in 1:ntabold # test for condition (14)
                        help2 = 0 
                        if l!=j && l!=k
                            for m in 1:help1 
                                if abs(Tab[l,C[m]])>EPS
                                    break 
                                else help2 += 1 
                                end
                            end
                        end 
                        if l == ntabold-1 
                            break 
                        end 
                        if help1 != help2 
                            B1[comb] = j 
                            B2[comb] = k
                            #comb += 1? 
                        end 
                    end
                end
            end
        end
    end 
    # if number of rows of new matrix Tab is greater than ntab, 
    # ntab increased and algorithm starts again
    if comb + n_zero > ntab 
        free(B1) 
        free(B2)
        free(Tab) 
        free(Z) 
        ntab += nflux 
        breaktest = 1 
        break 
    end
    
    #construction of the subsequent matrix Tab 
    for j in 1: comb 
        help1 = B1[j] 
        help2 = B2[j] 
        sgn = abs(Tab[help2,nflux+i])>abs(Tab[help1,nflux+i]) ? abs(Tab[help2,nflux+i]) : abs(Tab[help1,nflux+i])
        for k in 1:ntotal 
            Z[j+n_zero,k] = Tab[help1,k]*abs(Tab[help2,nflux+i])/sgn + Tab[help2,k]*abs(Tab[help1,nflux+i])/sgn 
        end
    end
    
    for j in 1:comb + n_zero 
        for k in 1:ntotal 
            Tab[j,k] = Z[j,k] 
            if k>nflux+i && abs(Tab[j,k]>EPS)
                breaktest = 1 
            end
        end
        if breaktest==0 # final tableau obtained
            break 
        end
    end
end

while breaktest == 1 && i<nmetab 
    comb += n_zero 
    # cancelling all row vectors in Tab the last element of which is zero,
    # storing them in Z, 
    # normalising the remaining rows of Tab by dividing them by their last element
    for i in 1:comb 
        if abs(Tab[i,nflux-1]) < EPS 
            for j in 1:nflux 
                Z[n_zero,j] = Tab[i,j] 
                n_zero += 1 
                comb -= 1 
                for j in i:comb 
                    for k in 1:nflux 
                        Tab[j,k] = Tab[j+1,k]
                        i -= 1
                    end
                end
            end
        end
    end
end
free(B1) 
free(B2)
free(Tab) 
free(Z) 
free(C) 
