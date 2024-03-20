#returns the connected component which includes a given node
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
        G[neighbor] = delete!(G[neighbor],node)
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

function branchorbits(graph::Dict)
    seqs = collect(keys(graph))
    seqs = filter(x->x.preperiod==0,seqs)
    seqs = filter(x->!('*' in x.items),seqs)
    return Set([Set(orbit(seq).items) for seq in seqs])
end
