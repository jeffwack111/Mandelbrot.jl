include("OrientTrees.jl")
include("../spiders/Spiders.jl")
include("../parameters/DynamicRays.jl")

function embednodes(OHT::OrientedHubbardTree)
    OZ = labelonezero(OHT)
    anglelist = allanglesof(OZ,OHT)

    criticalorbit = orbit(OHT.zero) #This is a sequence orbit

    theta = angleof(first(criticalanglesof(OZ,OHT)))
    c = parameter(standardspider(theta),250)
    print(c)

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
    
    merge!(rays,dynamicrays(c,1//2,100,10,40))
    merge!(rays,dynamicrays(c,0//1,100,10,40))

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

function standardedges(OHT::OrientedHubbardTree,zvalues)
    edges = edgeset(OHT.adj)

    edgevectors = []
    for edge in edges
        start = first(filter(z -> z!= OHT.zero,edge))
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

function edgepath(graph,start, finish)
    vertexpath = nodepath(graph, start,finish)
    return [Set(E) for E in partition(vertexpath,2,1)]
end

function refinedtree(OHT,zvalues,steps)
    c = zvalues[shift(OHT.zero)]
    E = standardedges(OHT,zvalues)
    for ii in 1:steps
        E = refinetree(OHT,c,E)
    end
    return E
end

function refinetree(OHT,c,edgevectors)
    newedges = []
    for edge in keys(edgevectors)
        start = edgevectors[edge][1]
        finish = first(filter(x->x!=start,edge))
        #determine the image of this edge
        image = longpath(OHT,edgevectors,shift(start),shift(finish))

        refinedimage = [image[1]]
        for point in image[2:end]
            push!(refinedimage, (refinedimage[end]+point)/2)
            push!(refinedimage,point)
        end

        x = sqrt(-c + refinedimage[1])

        if abs2(x-edgevectors[edge][2][1]) < abs2(-x-edgevectors[edge][2][1])
            push!(newedges,Pair(edge,(start,path_sqrt(-c .+ refinedimage))))
        else
            push!(newedges,Pair(edge,(start,-path_sqrt(-c .+ refinedimage))))
        end
    end
    return Dict(newedges)
end

function plotedges!(scene,edgevectors)
    for edge in edgevectors
        line = edge[2][2]
        lines!(scene,real.(line)/2,imag(line)/2,color = "black")
    end
    return scene
end

function plotedges(edgevectors)
    scene = Scene(size=(1000,1000),aspect = 1)
    return plotedges!(scene,edgevectors)
end

function embedanim(AIA::AngledInternalAddress,frames)
    OHT = OrientedHubbardTree(AIA)
    (E,c) = standardedges(OHT)

    edgelist = [E]
    for ii in 1:frames
        E = refinetree(OHT,c,E)
        push!(edgelist,E)
    end

    scene = Scene(size=(1000,1000),aspect = 1)
    record(scene,"embedding.gif",1:frames,framerate = 3) do ii
        empty!(scene)
        plotedges!(scene,edgelist[ii])
    end
end

function embedanim(angle::Rational,frames)
    AIA = AngledInternalAddress(angle)
    return embedanim(AIA,frames)
end

function showtree!(scene,angle::Rational)
    E = refinedtree(angle,8)
    return plotedges!(scene,E)
end

function showtree(angle::Rational)
    scene = Scene(size=(500,500))
    return showtree!(scene, angle)
end
    


