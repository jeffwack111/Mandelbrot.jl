using CairoMakie
using GraphMakie
using Graphs

include("HubbardTrees.jl")

function showtree((A, F,markedpoints))
    fig = Figure()

    ax1 = Axis(fig[1,1])
    graphplot!(ax1,SimpleGraph(A),nlabels = string.(collect(1:size(A)[1])))
    hidedecorations!(ax1); 
    hidespines!(ax1); 
    ax1.aspect = DataAspect(); 

    #construct the adjacency matrix of the dynamics graph
    B = zeros(Int64,size(A))
    for (ii,image) in enumerate(F)
        B[ii,image] = 1
    end
    ax2 = Axis(fig[1,2])
    graphplot!(ax2, SimpleDiGraph(B),nlabels = string.(collect(1:size(A)[1])) )
    hidedecorations!(ax2); 
    hidespines!(ax2); 
    ax2.aspect = DataAspect(); 

    return fig
end

function showtree(angle::Rational)
    return showtree(hubbardtree(angle))
end

function showtree(int::Vector{Int})
    return showtree(hubbardtree(int))
end