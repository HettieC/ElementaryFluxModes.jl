"""
$(TYPEDSIGNATURES)

Differentiate the usage of EFMs with respect to changes in the model parameters using implicit
differentiation of the Lagrangian.
Calculating the proportion of each EFM used in the optimal solution can be turned into the
following LP:
    maximise     ∑xᵢ
    subject to   Dx = 1
where D is the cost matrix associated to each EFM in each enzyme pool.
We then differentiate Lagrangian of this LP to calculate the differential of x by parameters.

Arguments:
- 'EFMs': a vector of dictionaries of EFMs
- 'parameters': vector of the model parameters
- `rid_pid`: dict of reaction ids => parameter ids
- `parameter_values`: dict of parameter symbol => float value
- `rid_gcounts`: dict of reaction id => gene product counts
- `capacity`: enzyme capacity limit of the enzyme constrained model
- `gene_product_molar_masses`: dict of gene product gene_product_molar_masses
Output:
Matrix Aᵢⱼ of the the differential of EFM i with respect to parameter j: d(EFMᵢ)/d(pⱼ)
"""
function differentiate_efm(
    EFMs::Vector{Dict{String,Float64}},
    parameters::Vector{FastDifferentiation.Node},
    rid_pid::Dict{String,FastDifferentiation.Node},
    parameter_values::Dict{Symbol,Float64},
    rid_gcounts::Dict{String,Dict{String,Float64}},
    capacity::Vector{Tuple{String,Vector{String},Float64}},
    gene_product_molar_masses::Dict{String,Float64},
    optimizer,
)
    D_eval =
        float.(
            _cost_matrix(
                EFMs,
                rid_pid,
                rid_gcounts,
                capacity,
                gene_product_molar_masses;
                evaluate = true,
                parameter_values,
            )
        )
    n_vars = size(D_eval, 2)

    efm_opt = JuMP.Model(optimizer)
    JuMP.@variable(efm_opt, z[1:n_vars])
    JuMP.@constraint(efm_opt, eq, D_eval * z == [1; 1])


    JuMP.@objective(efm_opt, Max, sum(z))
    JuMP.optimize!(efm_opt)

    x = JuMP.value.(efm_opt[:z])
    ν = JuMP.dual.(efm_opt[:eq])
    D =
        FastDifferentiation.Node.(
            _cost_matrix(EFMs, rid_pid, rid_gcounts, capacity, gene_product_molar_masses)
        )

    # define L, the gradient of the Lagrangian
    L(x, ν, parameters) = [
        ones(n_vars) + D' * ν
        D * x - ones(n_vars)
    ]
    # differentiate L wrt x,ν, the variables
    dl_vars = [
        SparseArrays.spzeros(n_vars, n_vars) D_eval'
        D_eval SparseArrays.spzeros(n_vars, n_vars)
    ]

    # differentiate L wrt parameters
    dL_params(x, ν, parameters) =
        FastDifferentiation.jacobian(L(x, ν, parameters), parameters)
    # substitute parameter values:
    dL_params_eval =
        FastDifferentiation.make_function(dL_params(x, ν, parameters), parameters)
    param_vals = float.(collect(values(parameter_values)))

    dx = -Array(dl_vars) \ dL_params_eval(param_vals)

    # note: dx[[3,4],:] gives the derivatives of the dual variables ν
    return dx[[1, 2], :]
end

export differentiate_efm

"""
$(TYPEDSIGNATURES)

Calculate a matrix of the cost vectors of each EFM to each constraint.
Entry (i,j) gives the total cost in constraint i to produce one unit objective flux through EFM j.
Cost is calculated as ∑w(i)V(j)/kcat, where the variables are:
- 'w(i)': fraction of the ith enzyme pool that one mole of the enzyme uses up
- 'V(j)': flux through the reaction in EFM j
- 'kcat': turnover number of the enzyme.
"""
function _cost_matrix(
    EFMs::Vector{Dict{String,Float64}},
    rid_pid,
    rid_gcounts,
    capacity::Vector{Tuple{String,Vector{String},Float64}},
    gene_product_molar_masses::Dict{String,Float64};
    evaluate = false,
    parameter_values = nothing,
)
    D = Matrix(undef, length(capacity), length(EFMs))
    for (i, (enzyme_group, enzymes, enzyme_bound)) in enumerate(capacity)
        for (j, efm) in enumerate(EFMs)
            D[i, j] = 0
            for (rid, gcount) in rid_gcounts
                pid = rid_pid[rid]

                # if genes not in pool i, skip
                all(((g, c),) -> g ∉ enzymes, gcount) && continue

                # otherwise, add the cost of this enzyme to the ith pool from the jth efm
                if !evaluate
                    D[i, j] += sum([
                        efm[rid] * gene_product_molar_masses[g] * c / (enzyme_bound * pid) for (g, c) in gcount if g ∈ enzymes
                    ])
                else
                    D[i, j] += Float64(
                        sum([
                            efm[rid] * gene_product_molar_masses[g] * c /
                            (enzyme_bound * parameter_values[Symbol(pid)]) for
                            (g, c) in gcount if g ∈ enzymes
                        ]),
                    )
                end
            end
        end
    end
    return D
end

differentiate_ofm(
    EFMs::Vector{Dict{String,Float64}},
    parameters::Vector{FastDifferentiation.Node},
    rid_pid::Dict{String,FastDifferentiation.Node},
    parameter_values::Dict{Symbol,Float64},
    rid_gcounts::Dict{String,Dict{String,Float64}},
    capacity::Vector{Tuple{String,Vector{String},Float64}},
    gene_product_molar_masses::Dict{String,Float64},
    optimizer,
) = differentiate_efm(
    EFMs,
    parameters,
    rid_pid,
    parameter_values,
    rid_gcounts,
    capacity,
    gene_product_molar_masses,
    optimizer,
)
