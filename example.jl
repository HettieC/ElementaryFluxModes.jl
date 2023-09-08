using ElementaryModes, Tulip
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

NS, pivotrows = rational_nullspace(S)
E, IFF, FAF = preprocessing(NS, pivotrows, Tulip)
E = main_programme(E, IFF, FAF, NS, S, Tulip)
