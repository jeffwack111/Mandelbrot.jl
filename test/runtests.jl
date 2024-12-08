using Spiders
using Test

@testset "Spiders.jl" begin
    @test KneadingSequence(3//5) == KneadingSequence(2//5)
end
