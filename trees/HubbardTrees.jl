#Algorithms from https://eudml.org/doc/283172
#Existence of quadratic Hubbard trees
#Henk Bruin, Alexandra Kafll, Dierk Schleicher

include("../sequences/AngleDoubling.jl") 

function hubbardtree(K::Sequence)
    starK = Sequence(pushfirst!(copy(K.items),'*'),K.preperiod+1)

    #We begin with the critical orbit
    markedpoints = copy(orbit(starK).items)

    #These points are added for orientation purposes
    #=
    push!(markedpoints,Sequence(['B'],0))
    push!(markedpoints,Sequence(['A','B'],1))
    push!(markedpoints,Sequence(['B','A'],1))
    push!(markedpoints,Sequence(['A'],0))
    =#

    n = length(markedpoints)

    results = fill(0,(n,n,n))

    runagain = true
    while runagain == true
        runagain = false
        #All triples of points on the critical orbit are run through the triod map
        for triple in subsets(collect(enumerate(markedpoints)),3)
            if results[triple[1][1],triple[2][1],triple[3][1]] == 0 #skipping the ones we've done already
                type,seq = iteratetriod(K,(triple[1][2],triple[2][2],triple[3][2]))

                if type == "branched"
                    if isempty(findall(x->x==seq,markedpoints))#If the branch point has not been found already,
                        push!(markedpoints,seq)                #add it to the list of marked points.
                        results = cat(results,1,dims=(1,2,3))  #We then must extend the results matrix to accomodate the new marked point
                        index = lastindex(markedpoints) 
                        runagain = true       
                    else
                        index = findall(x->x==seq,markedpoints)[1]
                    end
                elseif type == "flat"
                    index = triple[seq][1]
                end

                results[triple[1][1],triple[2][1],triple[3][1]] = index
            end
        end
    end

    #We've now found all the marked points of the tree, so d is the number of vertices
    d = length(markedpoints)

    M = fill([0],(d,d))

    for (ii,jj) in subsets(collect(1:d),2)
        middles = Int[]
        for kk in 1:d
            if kk != ii && kk != jj 
                triple = sort([ii,jj,kk])
                push!(middles,results[triple[1],triple[2],triple[3]])
            end
        end
        M[ii,jj] = middles
    end
   
    E = [Set{Sequence}() for ii in 1:d]

    for (ii,jj) in subsets(collect(1:d),2)
        pointsbetween = 0
        for mid in M[ii,jj]
            if mid != ii && mid !=jj
                pointsbetween += 1
            end
        end
        if pointsbetween == 0
            
            push!(E[ii],markedpoints[jj])
            push!(E[jj],markedpoints[ii])

        end
    end
    
    d = Dict(zip(markedpoints,E))
    return d            
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
            push!(chopList,0)

        elseif Set([triod[1].items[1],triod[2].items[1],triod[3].items[1]]) == Set(['A','B','*']) 
            if triod[1].items[1] == '*'
                indexofmiddle = 1
            elseif triod[2].items[1] == '*'
                indexofmiddle = 2
            elseif triod[3].items[1] == '*'
                indexofmiddle = 3
            end
            
            return "flat",indexofmiddle

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
        return "flat", 1
    elseif !(2 in chopList)
        return "flat", 2
    elseif !(3 in chopList)
        return "flat", 3
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

#returns the connected component which includes a given node
function component(graph::Dict,start)
    C = Set()
    activenodes = [start]
    while !isempty(activenodes)
        newactivenodes = []
        for node in activenodes
            for neighbor in graph[node]
                if !(neighbor in C)
                    push!(C,neighbor)
                    push!(newactivenodes,neighbor)
                end
            end
        end
        activenodes = newactivenodes
    end
    return C
end

function removenode(graph::Dict,node)
    G = copy(graph)
    neighbors = G[node]
    G = delete!(G,node)

    for neighbor in neighbors
        G[neighbor] = delete!(G[neighbor],node)
    end

    return G
end

function globalarms(graph::Dict,removed)
    neighbors = graph[removed]

    cutgraph = removenode(graph,removed)

    C = []

    for neighbor in neighbors
        push!(C,component(cutgraph,neighbor))
    end

    return C

end

function branchorbits(graph::Dict)
    seqs = collect(keys(graph))
    seqs = filter(x->x.preperiod==0,seqs)
    seqs = filter(x->!('*' in x.items),seqs)
    return Set([Set(orbit(seq).items) for seq in seqs])
end

#returns a set of points whose preimages form the entire tree. 
#This means every cycle will have one representative,
#and the representative from cycles of branch points will be the characteristic point.
function characteristicset(tree::Dict)
    orbs = branchorbits(tree)
    C = Sequence[]
    for orb in orbs
        #All the other points on this orbit and zero are in one arm, and the critical point is in another
        #4.1 page 34 TreesBook
    end
end



