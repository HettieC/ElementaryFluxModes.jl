# # Toy Model

using ElementaryFluxModes

using LinearAlgebra

# ## Load a simple model

# The code used to construct the model is located in `test/simple_model.jl`, but
# it is not shown here for brevity.

include("../../test/simple_model.jl"); #hide


model

# ## Prepare matrices for calculating EFMs

# Get the stoichiometric matrix of the model, this is what we use to find EFMs

N = AbstractFBCModels.stoichiometry(model)

# Take a rational nullspace of `N`

K = rational_nullspace(Matrix(N))[1]
K = round.(K, digits = 14)
# Permute the rows of K, so that it is in the form [I;K*]

order = Int64[]
rows_done = 1
for (i, row) in enumerate(eachrow(K))
    global rows_done
    rows_done > size(K, 2) && break
    if row == Matrix(I(size(K, 2)))[rows_done, :]
        push!(order, i)
        rows_done += 1
    end
end
append!(order, [i for i = 1:size(K, 1) if i ∉ order])
K = K[order, :]

# N also needs to be in same reaction order as K

N = N[:, order]

# ## Calculate the EFMs 

# Run the double description algorith on `N` and `K`

R = DDBinary(N, K)

# DDBinary returns a boolean matrix, flux values for
# each EFM then need to be calculated in E:

E = Matrix(undef, size(R, 1), size(R, 2))
for (i, r) in enumerate(eachcol(R))
    non_zero = findall(x -> x != 0, r)
    flux_ns = rational_nullspace(Matrix(N[:, non_zero]); tol = 1e-14)[1]
    mode = zeros(size(R, 1))
    for (j, x) in zip(non_zero, flux_ns)
        mode[j] = abs(x) < 1e-14 ? 0 : x
    end
    E[:, i] = mode
end

# Return the entries to the original reaction order

E = E[invperm(order), :]
