using ElementaryFluxModes
using Test

@testset "ElementaryFluxModes.jl" begin
    # standard
    @testset "Standard double dispatch EFMs" begin
        include("ddstandard.jl")
    end
    # binary trees...
end
