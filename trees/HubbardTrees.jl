#Algorithms from https://eudml.org/doc/283172
#Existence of quadratic Hubbard trees
#Henk Bruin, Alexandra Kafll, Dierk Schleicher

include("../sequences/AngleDoubling.jl") 
include("Graphs.jl")

function hubbardtree(K::Sequence)
    starK = Sequence(pushfirst!(copy(K.items),'*'),K.preperiod+1)
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
    return hubbardtree(Sequence(seq,0))
end

function iteratetriod(K::Sequence,triod::Tuple{Sequence,Sequence,Sequence})
    triodList = []
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
        return "branched",majorityvote(Sequence(triodList,preperiod))
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
    return Sequence(newitems,S.preperiod) 
end

function forwardimages(htree::Dict)
    return Dict([Pair(key,shift(key)) for key in keys(htree)])
end

#returns a set of points whose preimages form the entire tree. 
#This means every cycle will have one representative,
#and the representative from cycles of branch points will be the characteristic point.
function characteristicset(tree::Dict)
    orbs = orbits(tree)
    for orb in orbs
        #All the other points on this orbit and zero are in one arm, and the critical point is in another
        #4.1 page 34 TreesBook
    end
end



