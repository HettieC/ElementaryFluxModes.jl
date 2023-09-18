# tests for standard double dispatch

@testset "Toy model" begin
    const S = [
        1 -1 0 -1 0 0 0 0 0 0 0 0 0 0 ;
        0 1 -1 0 -1 1 0 0 0 0 0 0 0 0 ;
        0 0 0 1 1 -1 -1 -1 0 0 0 0 0 0 ;
        0 0 0 0 0 0 0 1 -1 1 -1 0 0 0 ;
        0 0 0 0 0 0 1 0 0 0 -1 0 0 1 ;
        0 0 0 0 0 0 0 0 0 0 1 -1 0 0 ;
        0 0 0 0 0 0 -1 0 0 0 0 1 -1 0 ;
        0 0 0 0 0 0 0 0 0 0 0 1 0 -1 ;
    ]
    @test true 
end

@testset "Core model" begin
    @test true
end

