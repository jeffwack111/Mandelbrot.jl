using CairoMakie
using NetworkLayout
using GraphMakie, Graphs

include("HubbardTrees.jl")

A = [1,3,6,7,8]

results, markedpoints = hubbardtree(A)

g = SimpleGraph(results)
layout = Stress()
f, ax, p = graphplot(g, layout=layout)
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f