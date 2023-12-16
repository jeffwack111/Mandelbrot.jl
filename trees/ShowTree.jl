using CairoMakie
using NetworkLayout
using GraphMakie, Graphs

include("HubbardTrees.jl")

function showtree((A, markedpoints))
    f, ax, = graphplot(SimpleGraph(A),nlabels = string.(collect(1:size(A)[1])))
    hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect(); 
    
    return f
end

function showtree(angle::Rational)
    return showtree(hubbardtree(angle))
end

function showtree(int::Vector{Int})
    return showtree(hubbardtree(int))
end