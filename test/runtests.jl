using ElementaryFluxModes
using Test

@testset "ElementaryFluxModes.jl" begin
    @testset "Standard double dispatch EFMs" begin
        include("ddstandard.jl")
    end
    @testset "Binary implementation of double dispatch EFMs" begin
        include("ddbinary.jl")
    end
    @testset "Toy model" begin
        include("../docs/src/1-toy-model.jl")
    end
end
