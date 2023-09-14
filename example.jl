using ElementaryModes, Tulip
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

R = DDStandard(S)
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


@time R = DDAlgorithm(Matrix(S))

findall(x -> x==false,E[:,2])



S = Matrix(deserialize("data/pgm_S"))
@time R = DDAlgorithm(S)

S = Matrix(deserialize("data/yeast_S"))
@time R = DDAlgorithm(S)
findall(x->x!=0,R)

round.(rational_nullspace(S)[1])