A = floor.(Int,[ 1.0   0.0  0.0  -1.0  0.0   0.0  0.0  1.0;
    -1.0   0.0  0.0   1.0  0.0   0.0  0.0  0.0  ; 
    -1.0   0.0  1.0   0.0  0.0   0.0  0.0  0.0  ;
    1.0   0.0  0.0   0.0  0.0   0.0  0.0  0.0  ;
    1.0   0.0  0.0   0.0  1.0  -1.0  0.0  1.0  
])
ns = nullspace(A)
ns = round.(nullspace(A),digits=12)

using Nemo
A = rationalize.(A)
r,ns = Nemo.nullspace_right_rational(ZZMatrix(A))

A*(Int.(Matrix(ns)))



#### test with ecoli core 
N = deserialize("data/pmodel_S_rat")
NZ = (sign.(N)) * (floor.(Int,abs.(N)))

