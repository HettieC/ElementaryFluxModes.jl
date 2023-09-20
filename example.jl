using ElementaryFluxModes
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
    0 0 0 0 0 0 0 0 0 0 0 1 0 -1.0 ;
]

ns = round.(rational_nullspace(S)[1])
nsrref = rref(ns')
R = nsrref'
E = DDStandard(S)

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

adjacency_test(E[:,1],E[:,2],E)

fixed_fluxes = [1,14]
flux_values = [2.0,2.0]
N1 = S[:,[i for i in 1:size(S,2) if i ∉ fixed_fluxes]]
N2 = S[:,[i for i in 1:size(S,2) if i ∈ fixed_fluxes]]
w = N2*flux_values
stoich = hcat(N1,w)

S_convex = fix_fluxes(S,fixed_fluxes,flux_values)
E = DDStandard(S_convex)


E = DDStandardFixedFluxes(S_convex, fixed_fluxes)

E_full = clean_DD_result(E, fixed_fluxes, flux_values)
v = E_full[[i for i in 1:size(E_full,1) if i ∉ fixed_fluxes],6]
N1*v == - N2*flux_values
N1*v

N1*v-N2*flux_values


#### test if rays are extreme 
A = hcat(N1, w)
zeta = [i for (i,x) in enumerate(E[:,8]) if x==0]
rank(vcat(A,I(size(A,2))[zeta,:])) == size(A,2) -1

rank(vcat(S,I(size(S,2))[zeta,:])) == size(S,2) -1


E = DDStandard(stoich)
E = clean_DD_result(E,[1,5],[2.0,2.0])

S = deserialize("data/EcoliCoreS")

@time R = DDStandard(Matrix(S))
findall(x -> x!=0,R)

findfirst(x->x!=0,E[2,:])

S = Matrix(deserialize("data/pgm_S"))
@time R = DDStandard(S)

S = Matrix(deserialize("data/yeast_S"))
@time R = DDStandard(Matrix(deserialize("data/yeast_S")))
findall(x->x!=0,R)

round.(rational_nullspace(S)[1])


S = deserialize("data/EcoliCore")
@time R = DDStandard(Matrix(S))

findall(x->x!=0,R)







function DDAlg(A::Matrix)
    d,n = size(A)
    Ad = A[1:d,1:d] # initial step
    ρ = Int64[]
    while ρ != collect(1:d)

    end

end