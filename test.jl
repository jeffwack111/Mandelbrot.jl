using Plots

include("Spiders.jl")

K = KneadingSequence(1//240)

L = standardlegs(K.orbit)
