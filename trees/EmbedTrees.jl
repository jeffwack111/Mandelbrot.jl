include("OrientTrees.jl")

function edgepath(graph,start, finish)
    vertexpath = nodepath(graph, start,finish)
    return [Set(E) for E in partition(vertexpath,2,1)]
end

function embednodes(OHT::OrientedHubbardTree)
    OZ = labelonezero(OHT)
    anglelist = allanglesof(OZ,OHT)

    criticalorbit = orbit(OHT.zero) #This is a sequence orbit

    theta = angleof(first(criticalanglesof(OZ,OHT)))
    c = spideriterate(theta,100)

    #critical orbit
    rays = Dict()
    #periodic branch points
    characteristicpoints = characteristicset(OHT)
    for point in characteristicpoints
        merge!(rays,dynamicrays(c,angleof(first(anglelist[point])),100,10,20))
    end

    #preperiod branch points. This is slow!
    for node in keys(OHT.adj)
        for angle in anglelist[node]
            if !(angle in keys(rays))
                theta = angleof(angle)
                push!(rays,Pair(theta,dynamicrays(c,theta,100,10,20)[theta]))
            end
        end
    end
    
    merge!(rays,dynamicrays(c,1//2,100,10,20))
    merge!(rays,dynamicrays(c,0//1,100,10,20))

    paramorbit = [0.0+0.0im]
    n = period(theta)
    for ii in 1:n-1
        push!(paramorbit,paramorbit[end]^2+c)
    end

    zvalues = Dict{Sequence,ComplexF64}()
    for node in keys(OHT.adj)
        if '*' in node.items #then it is in the critical orbit
            idx = findone(x->x=='*',node.items)
            push!(zvalues,Pair(node,paramorbit[mod1(n-idx+2,n)]))
        else
            list = [rays[angleof(angle)][end] for angle in anglelist[node]]
            push!(zvalues,Pair(node,sum(list)/length(list)))
        end
    end

    return zvalues
end

function standardedges(OHT::OrientedHubbardTree)
    zvalues = embednodes(OHT)
    edges = edgeset(OHT.adj)

    edgevectors = []
    for edge in edges
        start = first(edge)
        finish = first(filter(z -> z!= start,edge))
        push!(edgevectors,Pair(edge,(start,[zvalues[start],zvalues[finish]])))
    end
    return Dict(edgevectors)
end

function longpath(OHT, edgevectors, start, finish)
    nodes = nodepath(OHT.adj, start, finish)
    edges = edgepath(OHT.adj, start, finish)

    segments = [edgevectors[edge][2] for edge in edges]

    longpath = ComplexF64[]

    for (ii,edge) in enumerate(edges)
        if edgevectors[edge][1] == nodes[ii]
            if ii != 1
                append!(longpath,edgevectors[edge][2][2:end])
            else
                append!(longpath,edgevectors[edge][2])
            end
        else
            if ii != 1
                append!(longpath,reverse(edgevectors[edge][2])[2:end])
            else
                append!(longpath,reverse(edgevectors[edge][2]))
            end
        end
    end

    return longpath

end