using CairoMakie
using NetworkLayout
using GraphMakie, Graphs

include("HubbardTrees.jl")

function showtree(angle::Rational)

    A, markedpoints = hubbardtree(angle)

    g = SimpleGraph(A)
    layout = Stress()
    f, ax, p = graphplot(g, layout=layout)
    text!(stress(g);text = string.(collect(1:size(g)[1]).-1))
    hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f

    return f
end