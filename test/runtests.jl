using ElementaryFluxModes
using Test

import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
using JSONFBCModels
import FastDifferentiation as F
const Ex = F.Node
import DifferentiableMetabolism as D
using COBREXA
import COBREXA as X
using JSON
import Tulip as T
using CairoMakie

@testset "ElementaryFluxModes.jl" begin
    @testset "Binary implementation of double dispatch EFMs" begin
        include("ddbinary.jl")
    end
    @testset "Toy model EFMs" begin
        include("../docs/src/1-toy-model.jl")
    end
    @testset "Toy model differentiate" begin
        include("../docs/src/2-differentiate.jl")
    end
    @testset "E. coli core OFMs" begin
        include("../docs/src/3-ecoli-core.jl")
    end

end
