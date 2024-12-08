using Mandelbrot
using Test

@testset "Mandelbrot.jl" begin
    @test KneadingSequence(3//5) == KneadingSequence(2//5)
end
