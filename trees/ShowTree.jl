using CairoMakie
using NetworkLayout
using GraphMakie, Graphs

include("HubbardTrees.jl")

function showtree((A, markedpoints))

    g = SimpleGraph(A)
    layout = Stress()
    f, ax, p = graphplot(g, layout=layout)
    text!(stress(g);text = string.(collect(1:size(g)[1]).-1))
    hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); f

    return f
end

function showtree(angle::Rational)
    return showtree(hubbardtree(angle))
end

function showtree(int::Vector{Int})
    return showtree(hubbardtree(int))
end