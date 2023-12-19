using CairoMakie
using GraphMakie
using Graphs

include("HubbardTrees.jl")

function showtree((A, F,markedpoints))
    fig = Figure()

    n = size(A)[1]

    B = boundary(markedpoints)

    colors = fill("black",n)

    for index in B
        colors[index] = "red"
    end

    ax1 = Axis(fig[1,1])
    graphplot!(ax1,SimpleGraph(A),nlabels = string.(collect(1:n)),node_color = colors)
    hidedecorations!(ax1); 
    hidespines!(ax1); 
    ax1.aspect = DataAspect(); 

    #construct the adjacency matrix of the dynamics graph
    B = zeros(Int64,(n,n))
    for (ii,image) in enumerate(F)
        B[ii,image] = 1
    end
    ax2 = Axis(fig[1,2])
    graphplot!(ax2, SimpleDiGraph(B),nlabels = string.(collect(1:n)) )
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

function showbetatree((A, F,markedpoints))
    fig = Figure()

    ax1 = Axis(fig[1,1])
    graphplot!(ax1,SimpleGraph(A),nlabels = string.(collect(1:size(A)[1])),node_color = "red",edge_color = "red")
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

function showbetatree(angle::Rational)
    return showbetatree(betahubbardtree(angle))
end

function showbetatree(int::Vector{Int})
    return showbetatree(betahubbardtree(int))
end

function rootedtree((A, F,markedpoints))
    G = SimpleGraph(A)

    root = findall(x->x==1, F)[1]

    n = length(G.fadjlist)

    T = [[root],neighbors(G,root)]

    nadded = 1 + length(T[2])

    while nadded < n
        parents = T[end]
        children = []
        for vertex in parents
            for u in neighbors(G,vertex)
                if !(u in T[end-1])
                    push!(children,u)
                end
            end
        end
        push!(T,children)
        nadded += length(children)
    end

    Ngens = length(T)

    X = zeros(Float64,n)
    Y = zeros(Float64,n)
    
    for (gen,vertices) in enumerate(T)
        k = length(vertices)
        for (ii,u) in enumerate(vertices)
            X[u] = ii/(k+1)
            Y[u] = 1 - gen/(Ngens+1)
        end
    end

    B = boundary(markedpoints)

    colors = fill("black",n)

    for index in B
        colors[index] = "red"
    end


    f = Figure()

    ax1 = Axis(f[1, 1],aspect = 1)
    hidedecorations!(ax1)
    hidespines!(ax1)
    graphplot!(ax1,G,layout = Point.(X,Y),nlabels = string.(collect(1:n)),node_color = colors)

    B = zeros(Int64,(n,n))
    for (ii,image) in enumerate(F)
        B[ii,image] = 1
    end
    ax2 = Axis(f[1,2])
    graphplot!(ax2, SimpleDiGraph(B),nlabels = string.(collect(1:n)) )
    hidedecorations!(ax2); 
    hidespines!(ax2); 
    ax2.aspect = DataAspect(); 


    return f
end




    
