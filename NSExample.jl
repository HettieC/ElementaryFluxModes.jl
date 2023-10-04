A = [ 1.0   0.0  0.0  -1.0  0.0   0.0  0.0  1.0;
    -1.0   0.0  0.0   1.0  0.0   0.0  0.0  0.0  ; 
    -1.0   0.0  1.0   0.0  0.0   0.0  0.0  0.0  ;
    1.0   0.0  0.0   0.0  0.0   0.0  0.0  0.0  ;
    1.0   0.0  0.0   0.0  1.0  -1.0  0.0  1.0  
 ]
ns = nullspace(A)
ns = round.(nullspace(A),digits=12)

using Nemo
ns = Nemo.nullspace_right_rational(ZZMatrix(A))

## step 1 make matrix rational 
