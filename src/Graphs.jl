#implement adj and you should be able to use all these functions?
abstract type Graph end

function Base.getindex(G::Graph,idx)
    return G.adj[idx]
end

function vertices(G::Graph)
    return Set(keys(G.adj))
end

Base.show(io::IO,G::Graph) = display(G.adj)

function adjlist(graph)
    #first, assign an index to each point
    nodes = collect(keys(graph))
    E = Vector{Int}[]
    for node in nodes
        g = []
        for neighbor in graph[node]
            push!(g,findone(x->x==neighbor, nodes))
        end
        push!(E,g)
    end
    return (E,nodes)
end

function component(graph::Dict,start)
    C = Set{Sequence}()
    activenodes = [start]
    while !isempty(activenodes)
        newactivenodes = []
        for node in activenodes
            for neighbor in graph[node]
                if !(neighbor in C)
                    push!(C,neighbor)
                    push!(newactivenodes,neighbor)
                end
            end
        end
        activenodes = newactivenodes
    end
    return Pair(start,C)
end

function addto(graph,node,new)
    G = deepcopy(graph)
    push!(G[node],new)
    push!(G,Pair(new,Set([node])))
    return G
end

function addbetween(graph,A,B,new)
    G = deepcopy(graph)
    if !(B in G[A])
        error("A and B are not adjacent")
    end
    push!(G,Pair(new,Set([A,B])))
    push!(G[A],new)
    push!(G[B],new)
    delete!(G[A],B)
    delete!(G[B],A)
    return G
end

function removenode(graph::Dict,node)
    G = deepcopy(graph)
    neighbors = G[node]
    G = delete!(G,node)

    for neighbor in neighbors
        filter!(x-> x!= node,G[neighbor])
    end 

    return G
end

function globalarms(graph::Dict,removed)
    neighbors = graph[removed]

    cutgraph = removenode(graph,removed)

    C = []

    for neighbor in neighbors
        push!(C,component(cutgraph,neighbor))
    end

    return Dict(C)

end

function neighbortowards(graph,base,target)
    glarms = globalarms(graph,base)
    for (key,value) in glarms
        if key == target || target in value
            return key
        end
    end
end

function edgeset(graph::Dict)
    return Set(vcat([[Set([key,value]) for value in graph[key]] for key in keys(graph)]...))
end


