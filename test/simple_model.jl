
import AbstractFBCModels
const AC = AbstractFBCModels.CanonicalModel

rs = [
    AC.Reaction(;
        name = "r1",
        lower_bound = 0,
        upper_bound = 1000,
        stoichiometry = Dict("m1" => 1.0),
    )
    AC.Reaction(;
        name = "r2",
        lower_bound = 0,
        upper_bound = 1000,
        stoichiometry = Dict("m2" => 1),
    )
    AC.Reaction(;
        name = "ATPM",
        lower_bound = 0,
        upper_bound = 1000,
        stoichiometry = Dict("m4" => -1),
    )
    AC.Reaction(;
        name = "r3",
        lower_bound = 0,
        upper_bound = 1000,
        stoichiometry = Dict("m1" => -1, "m2" => -1, "m3" => 1),
        gene_association_dnf = [["g1"]],
    )
    AC.Reaction(;
        name = "r4",
        lower_bound = 0,
        upper_bound = 1000,
        stoichiometry = Dict("m3" => -1, "m4" => 1),
        gene_association_dnf = [["g2"]],
    )
    AC.Reaction(;
        name = "r5",
        lower_bound = 0,
        upper_bound = 1000,
        stoichiometry = Dict("m2" => -1, "m4" => 1),
        gene_association_dnf = [["g3", "g4"]],
    )
    AC.Reaction(;
        name = "r6",
        lower_bound = 0,
        upper_bound = 1000,
        stoichiometry = Dict("m4" => -1),
        objective_coefficient = 1.0,
    )
]

ms = [
    AC.Metabolite(; name = "m1")
    AC.Metabolite(; name = "m2")
    AC.Metabolite(; name = "m3")
    AC.Metabolite(; name = "m4")
]

gs = [
    AC.Gene(; name = "g1")
    AC.Gene(; name = "g2")
    AC.Gene(; name = "g3")
    AC.Gene(; name = "g4")
]

model = AC.Model(
    reactions = Dict(x.name => x for x in rs),
    metabolites = Dict(x.name => x for x in ms),
    genes = Dict(x.name => x for x in gs),
)

gene_product_molar_masses = Dict("g1" => 1.0, "g2" => 1.0, "g3" => 1.0, "g4" => 1.0)

rid_kcat = Dict("r3" => 10.0, "r4" => 10.0, "r5" => 5.0)

capacity = [("pool_1", ["g1", "g2"], 10.0), ("pool_2", ["g3", "g4"], 10.0)]

rid_gcounts = Dict(
    "r3" => Dict("g1" => 1.0),
    "r4" => Dict("g2" => 1.0),
    "r5" => Dict("g3" => 1.0, "g4" => 1.0),
)
