using ElementaryFluxModes, Tulip
using LinearAlgebra, RowEchelon, SparseArrays
using Serialization

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

#### make irreversible


S = [
    1 -1 0 -1 0 0 0 0  0 0 0 0 ;
    0 1 -1 0 -1 0 0 0  0 0 0 0 ;
    0 0 0 1 1 -1 -1 0  0 0 0 0 ;
    0 0 0 0 0 0 1 -1 -1 0 0 0 ;
    0 0 0 0 0 1 0 0 -1 0 0 1 ;
    0 0 0 0 0 0 0 0 1 -1 0 0 ;
    0 0 0 0 0 -1 0 0 0 1 -1 0 ;
    0 0 0 0 0 0 0 0 0 1 0 -1.0 ;
]

reversible = [5,8]

S_irrev = make_all_irreversible(S,reversible)


E = DDStandard(S_irrev)
E = reversible_EFMs(E,reversible)

E[:,[i for i in size(E,2) if any(x->x!=0,E[:,i])]]


@time R = DDStandard(S)
E = checkR(R,S)


r1 = R[:,5]
r2 = R[:,10]
zeta = [i for i in 1:length(r1) if r1[i] == 0 && r2[i]==0]
rank(R[zeta,:]) == rank(R)-2


## test if r is an extreme ray 
r = R[:,7]
zeta = [i for (i,x) in enumerate(r) if x==0]

rank(R[zeta_bar,:]) == length(zeta_bar) -2



n - length(zeta) == length(zeta_bar)



S = deserialize("data/EcoliCoreS")

@time R = DDStandard(Matrix(S))
findall(x -> x!=0,R)



S = Matrix(deserialize("data/pgm_S"))
@time R = DDStandard(S)

S = Matrix(deserialize("data/yeast_S"))
@time R = DDStandard(Matrix(deserialize("data/yeast_S")))
findall(x->x!=0,R)

round.(rational_nullspace(S)[1])


S = deserialize("data/EcoliCore")
@time R = DDStandard(Matrix(S))

findall(x->x!=0,R)