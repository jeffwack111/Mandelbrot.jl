include("OrientTrees.jl")
include("../spidermap/SpiderMap.jl")
include("../parameters/DynamicRays.jl")

struct EmbeddedHubbardTree
    zero::Sequence
    adj::Dict{Sequence,Vector{Sequence}}
    boundary::Vector{Sequence}
    onezero::Dict{Sequence{Char}, Char}
    angle::Rational
    rays::Dict{Rational,Vector{ComplexF64}}
    parameter::ComplexF64
    vertices::Dict{Sequence, ComplexF64} #maps each vertex in the tree to a point in the plane
    edges::Dict{Set{Sequence{Char}}, Tuple{Sequence{Char}, Vector{ComplexF64}}} #maps each edge to an oriented polyline
end

function EmbeddedHubbardTree(OHT::OrientedHubbardTree)

    OZ = labelonezero(OHT)
    anglelist = allanglesof(OZ,OHT)

    criticalorbit = orbit(OHT.zero) #This is a sequence orbit

    theta = angleof(first(criticalanglesof(OZ,OHT)))
    c = parameter(standardspider(theta),250)
    println(c)

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
                phi = angleof(angle)
                push!(rays,Pair(phi,dynamicrays(c,phi,100,10,20)[phi]))
            end
        end
    end
    
    merge!(rays,dynamicrays(c,1//2,100,10,40))
    merge!(rays,dynamicrays(c,0//1,100,10,40))

    paramorbit = [0.0+0.0im]
    n = period(theta)
    #println(n)
    for ii in 1:n-1
        push!(paramorbit,paramorbit[end]^2+c)
    end
    #println(paramorbit)

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
    E = standardedges(OHT.adj,OHT.zero,zvalues)
    return EmbeddedHubbardTree(OHT.zero,OHT.adj,OHT.boundary,OZ,theta,rays,c,zvalues,E)
end

function Base.show(io::IO,H::EmbeddedHubbardTree)
    return println(io,"Embedded Hubbard tree of "*repr(H.angle)*" with "*repr(length(keys(H.adj)))*" vertices")
end

function standardedges(adj,zero,zvalues)
    edges = edgeset(adj)

    edgevectors = []
    for edge in edges
        start = first(filter(z -> z!= zero,edge))
        finish = first(filter(z -> z!= start,edge))
        push!(edgevectors,Pair(edge,(start,[zvalues[start],zvalues[finish]])))
    end
    return Dict(edgevectors)
end

function longpath(EHT, start, finish)
    nodes = nodepath(EHT.adj, start, finish)
    edges = edgepath(EHT.adj, start, finish)

    edgevectors = EHT.edges

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

function refinetree!(EHT::EmbeddedHubbardTree)

    newedges = []
    for edge in keys(EHT.edges)
        start = EHT.edges[edge][1]
        finish = first(filter(x->x!=start,edge))
        #determine the image of this edge
        image = longpath(EHT,shift(start),shift(finish))

        refinedimage = [image[1]]
        for point in image[2:end]
            push!(refinedimage, (refinedimage[end]+point)/2)
            push!(refinedimage,point)
        end

        c = EHT.parameter
        x = sqrt(-c + refinedimage[1])

        if abs2(x-EHT.edges[edge][2][1]) < abs2(-x-EHT.edges[edge][2][1])
            EHT.edges[edge] = (start,path_sqrt(-c .+ refinedimage))
        else
            EHT.edges[edge] = (start,-path_sqrt(-c .+ refinedimage))
        end
    end
    return EHT
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
    

function angleclusters(H::OrientedHubbardTree)
    charset = characteristicset(H)
    allangles = allanglesof(H)
    anglegroups = [allangles[branch] for branch in charset]
    return anglegroups
end