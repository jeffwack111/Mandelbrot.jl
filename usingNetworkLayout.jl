using CairoMakie
using NetworkLayout
using GraphMakie, Graphs

include("HubbardTrees.jl")

A = [1,2,4,5,7,13,14]

results, markedpoints = hubbardtree(A)

g = SimpleGraph(results)
layout = Stress()
f, ax, p = graphplot(g, layout=layout)
text!(stress(g);text = string.(collect(1:size(g)[1]).-1))
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f