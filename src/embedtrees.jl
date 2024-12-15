struct HyperbolicComponent
    criticalpoint::Sequence
    adj::Dict{KneadingSequence,Vector{KneadingSequence}}
    boundary::Vector{KneadingSequence}
    onezero::Dict{KneadingSequence, Union{Nothing,Digit{2}}}
    angle::Rational
    rays::Dict{BinaryExpansion,Vector{ComplexF64}}
    parameter::ComplexF64
    vertices::Dict{KneadingSequence, ComplexF64} #maps each vertex in the tree to a point in the plane
    edges::Dict{Set{KneadingSequence}, Tuple{KneadingSequence, Vector{ComplexF64}}} #maps each edge to an oriented polyline
end

function HyperbolicComponent(OHT::OrientedHubbardTree)

    OZ = labelonezero(OHT)
    anglelist = allanglesof(OZ,OHT)

    criticalorbit = orbit(OHT.criticalpoint) #This is a sequence orbit

    theta = Rational(first(criticalanglesof(OZ,OHT)))
    c = parameter(standardspider(theta),250)

    #critical orbit
    rays = Dict()
    #periodic branch points
    characteristicpoints = characteristicset(OHT)
    for point in characteristicpoints
        merge!(rays,dynamicrays(c,first(anglelist[point]),100,10,20))
    end
    
    #preperiod branch points. This is slow!
    for node in keys(OHT.adj)
        for angle in anglelist[node]
            if !(angle in keys(rays))
                push!(rays,Pair(angle,dynamicrays(c,angle,100,10,20)[angle]))
            end
        end
    end
    
    merge!(rays,dynamicrays(c,BinaryExpansion(1//2),100,10,40))
    merge!(rays,dynamicrays(c,BinaryExpansion(0//1),100,10,40))

    paramorbit = [0.0+0.0im]
    n = period(theta)
    #println(n)
    for ii in 1:n-1
        push!(paramorbit,paramorbit[end]^2+c)
    end
    #println(paramorbit)

    zvalues = Dict{Sequence,ComplexF64}()
    for node in keys(OHT.adj)
        if star() in node.items #then it is in the critical orbit
            idx = findone(x->x==star(),node.items)
            push!(zvalues,Pair(node,paramorbit[mod1(n-idx+2,n)]))
        else
            list = [rays[angle][end] for angle in anglelist[node]]
            push!(zvalues,Pair(node,sum(list)/length(list)))
        end
    end
    E = standardedges(OHT.adj,OHT.criticalpoint,zvalues)
    return HyperbolicComponent(OHT.criticalpoint,OHT.adj,OHT.boundary,OZ,theta,rays,c,zvalues,E)
end

function Base.show(io::IO,H::HyperbolicComponent)
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
    c = zvalues[shift(OHT.criticalpoint)]
    E = standardedges(OHT,zvalues)
    for ii in 1:steps
        E = refinetree(OHT,c,E)
    end
    return E
end

function refinetree!(EHT)

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

function labelonezero(OH::OrientedHubbardTree)
    htree = OH.adj
    zero = OH.criticalpoint
    K = shift(zero)
    
    #We now want to traverse the boundary, assigning 1s and 0s to the entire tree as we go.
    #We will give the symbols to the neighbors at each vertex, so that in place of the Vector{Sequence{Char}}, we have a Vector{Tuple{Sequence{Char},{Char}}}.
    #Heading to minus beta from beta, 0 is on the left and 1 is on the right.
    #As we head along, we have a boundary point behind us and a boundary point ahead of us. 
    #The neighbors which are counterclockwise of the point behind and clockwise of the point ahead are in region 0
    boundary = OH.boundary
    OZ = Dict{KneadingSequence, Union{Nothing,Digit{2}}}([Pair(node,nothing) for node in boundary])

    lastnode = nothing

    for (ii,activenode) in enumerate(boundary[1:end-1])
        #first we label all of the activenode's neighbors.
        #If any of them are not boundary points, we label that entire branch.
        #what is the next node? 
        #It is the unique boundary point which is a neighbor of the current node and is not the last node.

        nextnode = boundary[ii+1]
        idx = findone(x-> x==nextnode,htree[activenode])

        num_neighbors = length(htree[activenode])

        idx = mod1(idx+1,num_neighbors)
        sym = Digit{2}(1)
        while htree[activenode][idx] != nextnode
            if htree[activenode][idx] == lastnode
                sym = Digit{2}(0)
            else
                push!(OZ, Pair(htree[activenode][idx],sym))
                #go down the rabbit hole
                keyvalue = component(removenode(htree,activenode),htree[activenode][idx])
                for child in push!(keyvalue[2],keyvalue[1])
                    push!(OZ,Pair(child,sym))
                end
            end
            idx = mod1(idx+1,num_neighbors)
        end
        lastnode = activenode
        activenode = nextnode
    end
            
    return OZ
end

function anglesofbranch(OZ,OH,node)
    htree = OH.adj
    boundary = OH.boundary
    l = node.preperiod
    initialneighbors = collect(htree[node])
    if length(initialneighbors) < 3
        error("$node is not a branch point")
    end
    neighbors = copy(initialneighbors)
    angles = [Digit{2}[] for x in neighbors]
    visited = []
    while !((node,neighbors) in visited) || isempty(visited)
        push!(visited,(node,copy(neighbors)))
        if OZ[node] === nothing #We are on the boundary and need to use neighbor info
            for (ii, neighb) in enumerate(neighbors)
                if OZ[neighb] === nothing
                    #Are we before or after the node on the boundary?
                    nodeidx = findone(x->x==node,boundary)
                    neighbidx = findone(x->x==neighb,boundary)
                    if nodeidx > neighbidx
                        push!(angles[ii],Digit{2}(0))
                    else
                        push!(angles[ii],Digit{2}(1))
                    end
                else
                    push!(angles[ii],OZ[neighb])
                end
            end
        else #We are not on the boundary and need no neighbor info
            for angle in angles
                push!(angle,OZ[node])
            end
        end
        
        node = shift(node)
        for (ii, neighb) in enumerate(neighbors)
            neighbors[ii] = neighbortowards(htree,node,shift(neighb))
        end
    end
    thetas = [BinaryExpansion(digits,l) for digits in angles]
    return thetas
end

function criticalanglesof(OZ,OH)
    H = OH.adj
    zero = OH.criticalpoint
    node = shift(zero)
    boundary = OH.boundary

    initialneighb = first(collect(H[node]))
    neighb = initialneighb
    angles = [Digit{2}[] for x in 1:2]
    visited = []
    while !(node in visited)
        push!(visited,node)
        if OZ[node] === nothing #We are on the boundary and need to use neighbor info
            for ii in 1:2
                if OZ[neighb] === nothing
                    #Are we before or after the node on the boundary?
                    nodeidx = findone(x->x==node,boundary)
                    neighbidx = findone(x->x==neighb,boundary)
                    if nodeidx > neighbidx 
                        if ii == 1
                            push!(angles[ii],Digit{2}(0))
                        else
                            push!(angles[ii],Digit{2}(1))
                        end
                    else
                        if ii == 1
                            push!(angles[ii],Digit{2}(1))
                        else
                            push!(angles[ii],Digit{2}(0))
                        end
                    end
                else
                    push!(angles[ii],OZ[neighb])
                end
            end
        else #We are not on the boundary and need no neighbor info
            for angle in angles
                push!(angle,OZ[node])
            end
        end
        node = shift(node)
        neighb = neighbortowards(H,node,shift(neighb))
    end
    thetadigits = [BinaryExpansion(digits,node.preperiod) for digits in angles]
    return thetadigits
end

function criticalanglesof(aia::AngledInternalAddress)
    OH = OrientedHubbardTree(aia)
    OZ = labelonezero(OH)
    return criticalanglesof(OZ,OH)
end

function allanglesof(OZ,OH)
    #first we do the critical orbit
    pairs = Dict()

    H = OH.adj
    boundary = OH.boundary
    zero = OH.criticalpoint
    seq = shift(zero)

    angles = criticalanglesof(OZ,OH)
    while !(seq in keys(pairs))
        push!(pairs,Pair(seq,angles))
        seq = shift(seq)
        angles = [shift(angle) for angle in angles]
    end

    #next we do the periodic orbits
    for c in characteristicset(OH)
        angles = anglesofbranch(OZ,OH,c)
        seq = c
        while !(seq in keys(pairs))
            push!(pairs,Pair(seq,angles))
            seq = shift(seq)
            angles = [shift(angle) for angle in angles]
        end
    end

    push!(pairs,Pair(KneadingSequence(KneadingSymbol[KneadingSymbol('B')],0),[BinaryExpansion([Digit{2}(0)],0)]))
    push!(pairs,Pair(KneadingSequence(KneadingSymbol[KneadingSymbol('A'),KneadingSymbol('B')],1),[BinaryExpansion([Digit{2}(1),Digit{2}(0)],1)]))

    #last we do the preperiodic branch points. This is slow!
    for node in keys(H)
        if !(node in keys(pairs))
            push!(pairs,Pair(node,anglesofbranch(OZ,OH,node)))
        end
    end

    return pairs 
end

function allanglesof(theta::Rational)
    OH = OrientedHubbardTree(theta)
    OZ = labelonezero(OH)
    return allanglesof(OZ,OH)
end

function allanglesof(H::OrientedHubbardTree)
    OZ = labelonezero(H)
    return allanglesof(OZ,H)
end

function angle_echo(theta::Rational)
    aia = AngledInternalAddress(theta)
    return criticalanglesof(aia)
end

function valid_binary(theta)
    aia = AngledInternalAddress(theta)
    K = kneadingsequence(theta)
    H = orientedtree(aia)
    return binary(theta) in criticalanglesof(H,K)
end




