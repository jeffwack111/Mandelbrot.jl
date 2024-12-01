include("HubbardTrees.jl") 

const OrientedHubbardTree = OrientedGraph{KneadingSequence}

function forwardimages(htree)
    return Dict([Pair(key,[shift(key)]) for key in keys(htree)])
end 

function preimages(htree)
    fimages = forwardimages(htree)
    return  Dict([Pair(key,filter(x -> key in fimages[x],keys(htree))) for key in keys(htree)])
end

function forwardorbit(htree,start)
    return component(forwardimages(htree), start)
end

function backwardsorbit(htree, start)
    return component(preimages(htree), start)
end

function grandorbit(htree, start)
    dynamicadjacency = mergewith(append!,forwardimages(htree),preimages(htree))
    p = component(dynamicadjacency, start)
    return Set(push!(p.second,p.first))
end

function dynamiccomponents(htree)
    return Set([grandorbit(htree,x) for x in keys(htree)])
end

#Then there are a unique point z ∈ {zk}n k=1 and two different components of T \ {z}
#such that the critical value is contained in one component and 0 and all other points
#zk not= z are in the other one
function ischaracteristic(htree,periodicpoint)
    glarms = globalarms(htree.adj,periodicpoint)

    if length(collect(keys(glarms))) < 3
        return false
    end

    forbit = forwardorbit(htree.adj,periodicpoint)[2]

    zero = htree.zero

    won = shift(zero)

    zarm = neighbortowards(htree.adj,periodicpoint,zero)

    for point in forbit
        if point != periodicpoint && !(point in glarms[zarm])
            return false
        end
    end
    
    if won in glarms[zarm]
        return false
    end

    return true

end


#This checks too many points with ischaracteristic. 
#The orbits of the tree should be divided up first, 
#then one characteristic point is found for each orbit
function characteristicset(H)
    nodes = collect(keys(H.adj))
    periodicnodes = filter(x -> x.preperiod == 0 , nodes)
    C = Sequence{KneadingSymbol}[]
    for node in periodicnodes
        if ischaracteristic(H,node)
            push!(C,node)
        end
    end
    #before returning this vector we must put it in order from 0 to the critical point.
    #This will result later in the correct numerators being matched up with the correct characteristic points
    #There should be a single element of C which is between 0 and all the other elements
    zero = H.zero
    crit_value = shift(zero)
    path = nodepath(H.adj,zero,crit_value)
    D = []
    for node in path
        if node in C
            push!(D,node)
        end
    end

    return D
end

function characteristicorbits(H)
    P = preimages(H.adj)
    return [Pair(c, Dict([Pair(x,P[x]) for x in component(P,c)[2]])) for c in characteristicset(H)]
end

function OrientedHubbardTree(AIA::AngledInternalAddress)
    #Construct the unoriented Hubbard tree from the internal address ignoring the angles
    H = HubbardTree(InternalAddress(AIA.addr))

    #Find the boundary between the '0' region and the '1' region as the path between the beta fixed point and its preimage
    (H,boundary) = addboundary(H)

    #Filter the angles for those describing the characteristic points
    charangles = filter(x -> denominator(x) != 2, AIA.angles)

    #create empty oriented tree. 
    orientedH = Dict{KneadingSequence, Vector{KneadingSequence}}()

    charorbits = characteristicorbits(H)
    for (node, ang) in zip(charorbits,charangles)
        merge!(orientedH,orientpreimages(H.adj,node[2],orientcharacteristic(H, node[1],ang)))
    end

    #Deal with the critical orbit
    zero = H.zero
    for node in orbit(zero).items
        push!(orientedH,Pair(node,collect(H.adj[node])))
    end
    
    mbeta = KneadingSequence(KneadingSymbol[A(),B()],1)
    #deal with beta orbit
    for node in orbit(mbeta).items
        push!(orientedH,Pair(node,collect(H.adj[node])))
    end

    if orientedH[zero][1].items[1] == B()
        circshift!(orientedH[zero],1)
    end

    return OrientedHubbardTree(orientedH)

end

function OrientedHubbardTree(angle::Rational)
    return OrientedHubbardTree(AngledInternalAddress(angle))
end
 
function orientnode(H,source::KneadingSequence,target::Pair{KneadingSequence, Vector{KneadingSequence}})
    
    sourcenode = source
    sourceneighbors = H[source]

    targetnode = target[1]
    targetneighbors = target[2]

    targetarms = globalarms(H,targetnode)
    orientedtargetarms = [push!(targetarms[arm],arm) for arm in targetneighbors]
    indexpairs = []

    for node in sourceneighbors
        image = shift(node)
        #find which arm of the target the image lays in
        armindex = findone(x-> image in x,orientedtargetarms)
        push!(indexpairs,Pair(node,armindex))
    end

    oriented = [keyvalue[1] for keyvalue in sort(indexpairs,by = x -> x[2])]

    return Pair(sourcenode,oriented)

end

#SKETCHY SKETCHY
function orientcharacteristic(H,node,ang)
    num = numerator(ang)
    den = denominator(ang)
    glarms = globalarms(H.adj,node)
    zero = H.zero
    per = zero.period
    c = shift(zero)
    q = length(glarms)
    n = node.period

    if den != q
        error("the denominator of the characteristic angle does not match the order of the branch point")
    end
    
    #the characteristic point has arms towards 1-n,1,1+n,1+2n, where n is the period of the characteristic point
    #the numerator tells us the arm towards 1 is in position num
    #with 0-indexing
    #1-n, 0%q
    #1,num%q
    #1+n,2*num
    indexpairs = [Pair(neighbortowards(H.adj,node,shiftby(c,per+n*x)),mod1(x*num,q)) for x in -1:(q-2)]
    oriented = [x[1] for x in sort(indexpairs, by=z->z[2])]
    return Pair(node,oriented)
end

function orientpreimages(H,P,target::Pair{KneadingSequence, Vector{KneadingSequence}})
    keyvaluepairs = Dict(target)
    activetargets = [target[1]]
    while !isempty(activetargets)
        newactivetargets = []
        for targ in activetargets
            for preim in P[targ]
                newkeyvalue = orientnode(H,preim,Pair(targ,keyvaluepairs[targ]))
                if !(preim in keys(keyvaluepairs))
                    push!(keyvaluepairs,newkeyvalue)
                    push!(newactivetargets,preim)
                end
            end
        end
        activetargets = newactivetargets
    end
    return keyvaluepairs
end

function addboundary(htree::HubbardTree)
    beta = KneadingSequence(KneadingSymbol[B()],0)
    mbeta = KneadingSequence(KneadingSymbol[A(),B()],1)
    H = extend(htree,beta)
    H = extend(H,mbeta)
    
    boundary = nodepath(H.adj,beta,mbeta)

    return (H,boundary)
end

function labelonezero(OH::OrientedHubbardTree)
    htree = OH.adj
    zero = OH.zero
    K = shift(zero)
    
    #We now want to traverse the boundary, assigning 1s and 0s to the entire tree as we go.
    #We will give the symbols to the neighbors at each vertex, so that in place of the Vector{Sequence{Char}}, we have a Vector{Tuple{Sequence{Char},{Char}}}.
    #Heading to minus beta from beta, 0 is on the left and 1 is on the right.
    #As we head along, we have a boundary point behind us and a boundary point ahead of us. 
    #The neighbors which are counterclockwise of the point behind and clockwise of the point ahead are in region 0
    boundary = OH.boundary
    OZ = Dict([Pair(node,'*') for node in boundary])
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
        sym = '1'
        while htree[activenode][idx] != nextnode
            if htree[activenode][idx] == lastnode
                sym = '0'
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
    angles = [Char[] for x in neighbors]
    visited = []
    while !((node,neighbors) in visited) || isempty(visited)
        push!(visited,(node,copy(neighbors)))
        if OZ[node] == '*' #We are on the boundary and need to use neighbor info
            for (ii, neighb) in enumerate(neighbors)
                if OZ[neighb] == '*'
                    #Are we before or after the node on the boundary?
                    nodeidx = findone(x->x==node,boundary)
                    neighbidx = findone(x->x==neighb,boundary)
                    if nodeidx > neighbidx
                        push!(angles[ii],'0')
                    else
                        push!(angles[ii],'1')
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
    thetas = [Sequence{Char}(digits,l) for digits in angles]
    return thetas
end

function criticalanglesof(OZ,OH)
    H = OH.adj
    zero = OH.zero
    node = shift(zero)
    boundary = OH.boundary

    initialneighb = first(collect(H[node]))
    neighb = initialneighb
    angles = [Char[] for x in 1:2]
    visited = []
    while !(node in visited)
        push!(visited,node)
        if OZ[node] == '*' #We are on the boundary and need to use neighbor info
            for ii in 1:2
                if OZ[neighb] == '*'
                    #Are we before or after the node on the boundary?
                    nodeidx = findone(x->x==node,boundary)
                    neighbidx = findone(x->x==neighb,boundary)
                    if nodeidx > neighbidx 
                        if ii == 1
                            push!(angles[ii],'0')
                        else
                            push!(angles[ii],'1')
                        end
                    else
                        if ii == 1
                            push!(angles[ii],'1')
                        else
                            push!(angles[ii],'0')
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
    thetadigits = [Sequence{Char}(digits,node.preperiod) for digits in angles]
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
    zero = OH.zero
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

    push!(pairs,Pair(Sequence("B",0),[Sequence("0",0)]))
    push!(pairs,Pair(Sequence("AB",1),[Sequence("10",1)]))

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






