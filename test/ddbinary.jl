# tests for binary double dispatch

@testset "Toy model" begin
    SToy = [
        1 -1 0 -1 0 0 0 0 0 0 0 0 0 0
        0 1 -1 0 -1 1 0 0 0 0 0 0 0 0
        0 0 0 1 1 -1 -1 -1 0 0 0 0 0 0
        0 0 0 0 0 0 0 1 -1 1 -1 0 0 0
        0 0 0 0 0 0 1 0 0 0 -1 0 0 1
        0 0 0 0 0 0 0 0 0 0 1 -1 0 0
        0 0 0 0 0 0 -1 0 0 0 0 1 -1 0
        0 0 0 0 0 0 0 0 0 0 0 1 0 -1
    ]
    KToy = [
        1.0 0.0 0.0 0.0 0.0 0.0
        0.0 1.0 0.0 0.0 0.0 0.0
        0.0 0.0 1.0 0.0 0.0 0.0
        0.0 0.0 0.0 1.0 0.0 0.0
        0.0 0.0 0.0 0.0 1.0 0.0
        0.0 0.0 0.0 0.0 0.0 1.0
        1.0 -1.0 -0.0 0.0 0.0 0.0
        0.0 -1.0 1.0 1.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0
        1.0 0.0 -1.0 0.0 0.0 0.0
        1.0 0.0 -1.0 0.0 -1.0 1.0
        1.0 0.0 -1.0 0.0 -1.0 1.0
        1.0 0.0 -1.0 0.0 -1.0 1.0
        1.0 0.0 -1.0 0.0 -1.0 1.0
    ]
    RToy = [
        1 0 0 1 1 1 1 0 1
        0 0 0 1 1 0 0 0 1
        0 0 0 1 0 1 0 0 0
        0 1 0 0 1 0 0 0 1
        0 0 0 0 0 0 1 1 1
        0 0 1 0 0 0 0 1 0
        1 0 0 0 0 1 1 0 0
        0 1 0 0 0 1 0 0 0
        0 0 0 0 0 0 0 0 0
        1 0 0 0 1 0 1 0 1
        1 0 1 0 1 0 0 0 0
        1 0 1 0 1 0 0 0 0
        1 0 1 0 1 0 0 0 0
        1 0 1 0 1 0 0 0 0
    ]
    @test all(col -> col ∈ eachcol(RToy), eachcol(DDBinary(SToy, KToy)))
end

@testset "Core model" begin
    @test true
end
