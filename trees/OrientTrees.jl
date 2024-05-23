include("HubbardTrees.jl")

function forwardimages(htree::Dict)
    return Dict([Pair(key,[shift(key)]) for key in keys(htree)])
end 

function preimages(htree::Dict)
    fimages = forwardimages(htree)
    return  Dict([Pair(key,filter(x -> key in fimages[x],keys(htree))) for key in keys(htree)])
end

function forwardorbit(htree::Dict,start)
    return component(forwardimages(htree), start)
end

function backwardsorbit(htree::Dict, start)
    return component(preimages(htree), start)
end

function grandorbit(htree::Dict, start)
    dynamicadjacency = mergewith(append!,forwardimages(htree),preimages(htree))
    p = component(dynamicadjacency, start)
    return Set(push!(p.second,p.first))
end

function dynamiccomponents(htree::Dict)
    return Set([grandorbit(htree,x) for x in keys(htree)])
end

#Then there are a unique point z ∈ {zk}n k=1 and two different components of T \ {z}
#such that the critical value is contained in one component and 0 and all other points
#zk not= z are in the other one
function ischaracteristic(htree,periodicpoint)
    glarms = globalarms(htree,periodicpoint)

    if length(collect(keys(glarms))) < 3
        return false
    end

    forbit = forwardorbit(htree,periodicpoint)[2]

    zero = first(filter(x->x.items[1]=='*',keys(htree)))

    won = shift(zero)

    zarm = neighbortowards(htree,periodicpoint,zero)

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
function characteristicset(tree::Dict)
    nodes = collect(keys(tree))
    periodicnodes = filter(x -> x.preperiod == 0 , nodes)
    C = Sequence{Char}[]
    for node in periodicnodes
        if ischaracteristic(tree,node)
            push!(C,node)
        end
    end
    return C
end

function orientedtree(AIA::AngledInternalAddress)
    H = hubbardtree(AIA.addr)

    (H,boundary) = addboundary(H)

    charpoints = characteristicset(H)
    nums = Int[]
    for angle in AIA.angles
        if denominator(angle) > 2
            push!(nums,numerator(angle))
        end
    end

    orientedH = Dict{Sequence{Char}, Vector{Sequence{Char}}}()

    for (node, num) in zip(charpoints,nums)
        merge!(orientedH,orientpreimages(H,orientcharacteristic(H, node,num)))
    end

    #Deal with the critical orbit
    zero = first(filter(x->x.items[1]=='*',keys(H)))
    for node in orbit(zero).items
        push!(orientedH,Pair(node,collect(H[node])))
    end
    
    mbeta = Sequence("AB",1)
    #deal with beta orbit
    for node in orbit(mbeta).items
        push!(orientedH,Pair(node,collect(H[node])))
    end

    if orientedH[zero][1].items[1] == 'B'
        circshift!(orientedH[zero],1)
    end

    return orientedH,boundary

end

function orientedtree(angle::Rational)
    return orientedtree(AngledInternalAddress(angle))
end

function orientnode(H,source::Sequence{Char},target::Pair{Sequence{Char}, Vector{Sequence{Char}}})
    
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
function orientcharacteristic(H,node,num)
    glarms = globalarms(H,node)
    zero = first(filter(x->x.items[1]=='*',keys(H)))
    per = period(zero)
    c = shift(zero)
    q = length(glarms)
    n = period(node)
    
    #the characteristic point has arms towards 1-n,1,1+n,1+2n, where n is the period of the characteristic point
    #the numerator tells us the arm towards 1 is in position num
    #with 0-indexing
    #1-n, 0%q
    #1,num%q
    #1+n,2*num
    indexpairs = [Pair(neighbortowards(H,node,shiftby(c,per+n*x)),mod1(x*num,q)) for x in -1:(q-2)]
    oriented = [x[1] for x in sort(indexpairs, by=z->z[2])]
    return Pair(node,oriented)
end

function orientpreimages(H,target::Pair{Sequence{Char}, Vector{Sequence{Char}}})
    P = preimages(H)
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

function isbetween(htree::Dict,a,b,c)
    zero = first(filter(x->x.items[1]=='*',keys(htree)))
    K = shift(zero)
    (type,vertex) = iteratetriod(K,(a,b,c))
    if type == "flat" && vertex == a
        return true
    else
        return false
    end
end

function addboundary(htree::Dict)
    beta = Sequence("B",0)
    mbeta = Sequence("AB",1)
    H = addsequence(htree,beta)
    H = addsequence(H,mbeta)
    
    boundary = [beta]
    node = beta
    while node != mbeta
        node = first(filter(x->isbetween(H,x,mbeta,node),H[node]))
        push!(boundary,node)  
    end

    return (H,boundary)
end

function labelonezero(htree::Dict{Sequence{Char}, Vector{Sequence{Char}}},boundary)

    zero = first(filter(x->x.items[1]=='*',keys(htree)))
    K = shift(zero)
    
    #We now want to traverse the boundary, assigning 1s and 0s to the entire tree as we go.
    #We will give the symbols to the neighbors at each vertex, so that in place of the Vector{Sequence{Char}}, we have a Vector{Tuple{Sequence{Char},{Char}}}.
    #Heading to minus beta from beta, 0 is on the left and 1 is on the right.
    #As we head along, we have a boundary point behind us and a boundary point ahead of us. 
    #The neighbors which are counterclockwise of the point behind and clockwise of the point ahead are in region 0

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

function anglesof((htree,boundary),node)
    OZ = labelonezero(htree,boundary)
    return anglesof(OZ,(htree,boundary),node)
end

#Depreciated. Make a function that calculates the angles of all the branch points of a tree 
#Importantly, this function will already know the angle of the tree
function anglesof(OZ,(htree,boundary),node)
   
    if length(collect(htree[node])) > 1
        initialneighbors = collect(htree[node])
        neighbors = copy(initialneighbors)
        angles = [Char[] for x in neighbors]
        visited = []
        while !((node in visited) && neighbors == initialneighbors) || isempty(visited)
            push!(visited,node)
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
    else #We are on a leaf
        initialneighb = first(collect(htree[node]))
        neighb = initialneighb
        angles = [Char[] for x in 1:2]
        visited = []
        while !((node in visited) && neighb == initialneighb) || isempty(visited)
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
            neighb = neighbortowards(htree,node,shift(neighb))
        end
    end
    thetadigits = [Sequence{Char}(digits,node.preperiod) for digits in angles]
    thetas = sort([angleof(t) for t in thetadigits])
    return thetas
end

function anglesof(aia::AngledInternalAddress)
    H,boundary = orientedtree(aia)
    OZ = labelonezero(H,boundary)

    zero = first(filter(x->x.items[1]=='*',keys(H)))
    node = shift(zero)

    initialneighb = first(collect(H[node]))
    neighb = initialneighb
    angles = [Char[] for x in 1:2]
    visited = []
    while !((node in visited) && neighb == initialneighb) || isempty(visited)
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
    thetas = sort([angleof(t) for t in thetadigits])
    return thetas
end

function angle_echo(theta::Rational)
    aia = AngledInternalAddress(theta)
    return anglesof(aia)
end

function valid_binary(theta)
    aia = AngledInternalAddress(theta)
    K = kneadingsequence(theta)
    H = orientedtree(aia)
    return binary(theta) in anglesof(H,K)
end





