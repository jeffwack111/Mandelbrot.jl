#Algorithms from https://eudml.org/doc/283172
#Existence of quadratic Hubbard trees
#Henk Bruin, Alexandra Kafll, Dierk Schleicher

include("../sequences/AngleDoubling.jl") 
include("Graphs.jl")

function hubbardtree(K::Sequence{Char})
    starK = Sequence{Char}(pushfirst!(copy(K.items),'*'),K.preperiod+1)
    #We begin with the critical orbit
    markedpoints = copy(orbit(starK).items)

    H = Dict([Pair(starK,Set([K])),Pair(K,Set([starK]))])
    #println(H)
    for point in markedpoints[3:end]
        if !(point in keys(H)) 
            #println("adding $point")
            H = addsequence(H,K,(starK,deepcopy(H[starK])),point)
            #println("result is $H")
        end
    end
    return H
end

function addsequence(Htree,Kseq,(A,Bset),newpoint)
    if newpoint in collect(keys(Htree))
        print("warning: point already in tree")
        return Htree
    end
    B = pop!(Bset)
    
    #println("running triod($A,$B,$newpoint)")
    (result,point) = iteratetriod(Kseq,(A,B,newpoint))
    #println("$result, $point")
    if result == "branched"
        return addto(addbetween(Htree,A,B,point),point,newpoint)
    elseif result == "flat" && point == newpoint
        return addbetween(Htree,A,B,newpoint)
    else #then result is flat and the middle is A or B 
        if length(collect(Htree[point]))==1 #then point is a leaf
            #println("LEAF")
            return addto(Htree,point,newpoint)
        elseif point == B #if B is middle then B points down the right path
            forwardneighbors = delete!(deepcopy(Htree[B]),A)
            #println("STEP")
            return addsequence(Htree,Kseq,(B,forwardneighbors),newpoint)
        elseif isempty(Bset)
            #println("ATTATCH")
            return addto(Htree,point,newpoint)
        else #then we need to try a new path
            #println("ROTATE")
            return addsequence(Htree,Kseq,(A,Bset),newpoint)
        end
    end   
end

function hubbardtree(angle::Rational)
    return hubbardtree(kneadingsequence(angle))
end

function hubbardtree(internaladdress::Vector{Int})
    K = kneadingsequence(internaladdress)
    seq = copy(K.items)
    seq[end] = '*'
    return hubbardtree(Sequence{Char}(seq,0))
end

function iteratetriod(K::Sequence,triod::Tuple{Sequence,Sequence,Sequence})
    triodList = Tuple{Sequence,Sequence,Sequence}[]
    chopList = []

    while isempty(findall(x->x==triod,triodList))
        push!(triodList,triod)

        if triod[1].items[1] == triod[2].items[1] && triod[1].items[1] == triod[3].items[1]
            triod = (shift(triod[1]),shift(triod[2]),shift(triod[3]))

        elseif Set([triod[1].items[1],triod[2].items[1],triod[3].items[1]]) == Set(['A','B','*']) 
            if triod[1].items[1] == '*'
                middle = triodList[1][1]
            elseif triod[2].items[1] == '*'
                middle = triodList[1][2]
            elseif triod[3].items[1] == '*'
                middle = triodList[1][3]
            end
            
            return "flat",middle

        elseif triod[1].items[1] == triod[2].items[1] 
            triod = (shift(triod[1]),shift(triod[2]),K)
            push!(chopList,3)

        elseif triod[1].items[1] == triod[3].items[1]
            triod = (shift(triod[1]),K,shift(triod[3]))
            push!(chopList,2)

        elseif triod[2].items[1] == triod[3].items[1]
            triod = (K,shift(triod[2]),shift(triod[3]))
            push!(chopList,1)

        end

    end

    preperiod = findall(x->x==triod,triodList)[1] - 1

    if (1 in chopList) && (2 in chopList) && (3 in chopList) 
        return "branched",majorityvote(Sequence{Tuple{Sequence,Sequence,Sequence}}(triodList,preperiod))
    elseif !(1 in chopList)
        return "flat", triodList[1][1]
    elseif !(2 in chopList)
        return "flat", triodList[1][2]
    elseif !(3 in chopList)
        return "flat", triodList[1][3]
    end
end

function majorityvote(arms::Tuple{Sequence,Sequence,Sequence})
    if arms[1].items[1] == arms[2].items[1] 
        return arms[1].items[1]
    elseif arms[1].items[1] == arms[3].items[1]
        return arms[1].items[1]
    elseif arms[2].items[1] == arms[3].items[1]
        return arms[3].items[1]
    end
end

function majorityvote(S::Sequence)
    newitems = Char[]
    for triod in S.items
        push!(newitems,majorityvote(triod))
    end
    return Sequence{Char}(newitems,S.preperiod) 
end

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

function embed(AIA::AngledInternalAddress)
    H = hubbardtree(AIA.addr)
    charpoints = characteristicset(H)
    nums = Int[]
    for angle in AIA.angles
        if denominator(angle) > 2
            push!(nums,numerator(angle))
        end
    end

    orientedH = Dict()

    for (node, num) in zip(charpoints,nums)
        merge!(orientedH,orientpreimages(H,orientcharacteristic(H, node,num)))
    end

    zero = first(filter(x->x.items[1]=='*',keys(H)))
    for node in orbit(zero).items
        push!(orientedH,Pair(node,collect(H[node])))
    end

    if orientedH[zero][1].items[1] == 'B'
        circshift!(orientedH[zero],1)
    end

    return orientedH

end

function embed(angle::Rational)
    return embed(AngledInternalAddress(angle))
end


function orient(H,source::Sequence{Char},target::Pair{Sequence{Char}, Vector{Sequence{Char}}})
    
    sourcenode = source
    sourceneighbors = H[source]

    targetnode = target[1]
    targetneighbors = target[2]

    targetarms = globalarms(H,targetnode)
    orientedtargetarms = [push!(targetarms[arm],arm) for arm in targetneighbors]
    #println(orientedtargetarms)
    indexpairs = []

    for node in sourceneighbors
        
        image = shift(node)
        #println(image)
        #find which arm of the target the image lays in
        armindex = findone(x-> image in x,orientedtargetarms)
        push!(indexpairs,Pair(node,armindex))
    end

    oriented = [keyvalue[1] for keyvalue in sort(indexpairs,by = x -> x[2])]

    return Pair(sourcenode,oriented)

end

function orientcharacteristic(H,node,num)
    glarms = globalarms(H,node)
    zero = first(filter(x->x.items[1]=='*',keys(H)))
    c = shift(zero)
    degree = length(glarms)
    k = period(node)
    #the characteristic point has arms towards 0,1,1+k,1+2k, where k is the period of the characteristic point
    #the numerator tells us the arm towards 1 is in position num
    #with 0-indexing
    #0, 0%n
    #1,num%n
    #1+k,2*num%n
    indexpairs = [Pair(neighbortowards(H,node,shiftby(c,k*x)),mod((x+1)*num,degree)) for x in 0:degree-1]
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
                newkeyvalue = orient(H,preim,Pair(targ,keyvaluepairs[targ]))
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

function findone(f,A)
    list = findall(f,A)
    if length(list) == 0
        error("none found!")
    elseif length(list) > 1
        error("more than one found!")
    else
        return first(list)
    end
end