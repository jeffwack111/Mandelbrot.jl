#Algorithms from https://eudml.org/doc/283172
#Existence of quadratic Hubbard trees
#Henk Bruin, Alexandra Kafll, Dierk Schleicher

include("../sequences/AngleDoubling.jl") 
include("Graphs.jl")

abstract type AbstractHubbardTree 
end

struct HubbardTree <: AbstractHubbardTree
    zero::Sequence #Rename this to critical? Remove it and use H["*"] everywhere?
    adj::Dict{Sequence,Set{Sequence}}
end

function Base.getindex(H::HubbardTree, str::String)
    list = []
    for (K,neighbors) in pairs(H.adj)
        if agrees(K,str)
            newneighbors = []
                for N in neighbors
                    if agrees(N,str)
                        push!(newneighbors,N)
                    end
                end
            push!(list,Pair(K,newneighbors))
        end
    end
    return Dict(list)
end

function agrees(K::Sequence,str::String)
    for (ii,item) in enumerate(str)
        if item !== K[ii] 
            return false
        end
    end
    return true
end

function HubbardTree(K::Sequence{Char})
    starK = prepend(K,'*')
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
    return HubbardTree(starK,H)
end

function HubbardTree(internaladdress::Vector{Int})
    K = kneadingsequence(internaladdress)
    seq = copy(K.items)
    seq[end] = '*'
    return HubbardTree(Sequence{Char}(seq,0))
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

function extend(H::HubbardTree,new::Sequence)
    Z = H.zero
    K = shift(Z)
    return HubbardTree(Z,addsequence(H.adj,K,(Z,deepcopy(H.adj[Z])),new))
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

#This function is a wrapper for findall which throws an error if there is not exactly one element satisfying f
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

function nodepath(graph, start, finish)
    path = [start]
    node = start
    while node != finish
        node = first(filter(x->isbetween(graph,x,finish,node),graph[node]))
        push!(path,node)  
    end
    return path
end
