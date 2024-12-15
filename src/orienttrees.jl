
"A combinatorial Hubbard tree with cyclic order at each vertex. 
It can be constructed from an AngledInternalAddress
or a RationalAngle, and it can be used to construct a HyperbolicComponent"
struct OrientedHubbardTree <: AbstractHubbardTree
    adj::Dict{Sequence,Vector{Sequence}}
    criticalpoint::KneadingSequence
    boundary::Vector{KneadingSequence}
end

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

function grandorbit(htree, start)
    dynamicadjacency = mergewith(append!,forwardimages(htree),preimages(htree))
    p = component(dynamicadjacency, start)
    return Set(push!(p.second,p.first))
end

#Should be used to plot tree dynamics?
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

    zero = htree.criticalpoint

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
    zero = H.criticalpoint
    crit_value = shift(zero)
    path = nodepath(H,zero,crit_value)
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
    zero = H.criticalpoint
    for node in orbit(zero).items
        push!(orientedH,Pair(node,collect(H.adj[node])))
    end
    
    mbeta = KneadingSequence(KneadingSymbol[KneadingSymbol('A'),KneadingSymbol('B')],1)
    #deal with beta orbit
    for node in orbit(mbeta).items
        push!(orientedH,Pair(node,collect(H.adj[node])))
    end

    if orientedH[zero][1].items[1] == KneadingSymbol('B')
        circshift!(orientedH[zero],1)
    end

    return OrientedHubbardTree(orientedH,zero,boundary)

end

function OrientedHubbardTree(angle::Rational)
    return OrientedHubbardTree(AngledInternalAddress(RationalAngle(angle)))
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
    zero = H.criticalpoint
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
    beta = KneadingSequence(KneadingSymbol[KneadingSymbol('B')],0)
    mbeta = KneadingSequence(KneadingSymbol[KneadingSymbol('A'),KneadingSymbol('B')],1)
    H = extend(htree,beta)
    H = extend(H,mbeta)
    
    boundary = nodepath(H,beta,mbeta)

    return (H,boundary)
end