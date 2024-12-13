using ElementaryFluxModes
using Test

@testset "ElementaryFluxModes.jl" begin
    @testset "Binary implementation of double dispatch EFMs" begin
        include("ddbinary.jl")
    end
    @testset "Toy model" begin
        include("../docs/src/1-toy-model.jl")
    end
    @testset "E coli" begin 
        include("../docs/src/2-ecoli-core.jl")
    end
    @testset "Differentiate EFMs" begin 
        include("../docs/src/3-differentiate-efms.jl")
    end
end
