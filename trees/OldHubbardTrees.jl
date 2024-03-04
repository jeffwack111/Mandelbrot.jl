#Algorithms from https://eudml.org/doc/283172
#Existence of quadratic Hubbard trees
#Henk Bruin, Alexandra Kafll, Dierk Schleicher

include("../sequences/AngleDoubling.jl") 

function hubbardtree(K::Sequence)
    starK = Sequence(pushfirst!(copy(K.items),'*'),K.preperiod+1)

    #We begin with the critical orbit
    markedpoints = circshift(copy(orbit(starK).items),-1)

    #These points are added for orientation purposes
    push!(markedpoints,Sequence(['B'],0))
    push!(markedpoints,Sequence(['A','B'],1))
    push!(markedpoints,Sequence(['B','A'],1))
    push!(markedpoints,Sequence(['A'],0))

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
   
    #E = [Set{Sequence}() for ii in 1:d]
    E = [Set{Int}() for ii in 1:d]

    for (ii,jj) in subsets(collect(1:d),2)
        pointsbetween = 0
        for mid in M[ii,jj]
            if mid != ii && mid !=jj
                pointsbetween += 1
            end
        end
        if pointsbetween == 0
            
            #push!(E[ii],markedpoints[jj])
            #push!(E[jj],markedpoints[ii])
            
            push!(E[ii],jj)
            push!(E[jj],ii)
        end
    end
    
    #We now calculate the vector describing the dynamics of the tree.
    #The nth entry of the vector hold the index of the image of the nth sequence under the shift map
    F = Int[]
    for seq in markedpoints
        append!(F,findall(x->x==shift(seq),markedpoints)[1])
    end

    return (E, F, markedpoints)
    

    #=
    d = Dict(zip(markedpoints,E))
    return d
    =#
                    
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

function boundary(markedpoints)
    beta = Sequence(['B'],0)
    mbeta = Sequence(['A','B'],1)

    B = Int[]

    for (ii,point) in enumerate(markedpoints)
        type,seq = iteratetriod(markedpoints[1],(beta,mbeta,point))

        if type == "flat"
            if seq == 3
                push!(B,ii)
            end
        end
    end

    return B
end

function embedtree(AIA::AngledInternalAddress)
    (E,F,markedpoints) = hubbardtree(AIA.addr)

    R = [collect(x) for x in E]

    numerators = Int[]
    for angle in AIA.angles
        if denominator(angle) != 2
            push!(numerators, numerator(angle))
        end
    end

    #First we want to use the numerators 
    #to assign cyclic order to characteristic points
    C = characteristicpoints((E,F,markedpoints))
    
    if length(C) != length(numerators)
        error("mismatch between #(numerators) and #(characteristicpoints)")
    end
    
    z = findall(x->x==1,F)[1]

    for (jj,chpoint) in enumerate(C)
        
        GLARM = globalarms(E, chpoint)

        d = length(E[chpoint])
        order = zeros(Int64,d)

        #what is the period of this critical point?
        k = period(markedpoints[chpoint])

        #which is the arm towards the critical point?
        for t in 1:d
            for arm in GLARM
                if (1+(t-1)*k) in arm
                    order[mod1(numerators[jj]*t,d)] = arm[2]
                end
            end
        end

        R[chpoint] = order  
        
        
        #now we have to go through all the preimages of this characteristic point

        activenodes = [chpoint]
        while !isempty(activenodes)
            newactivenodes = []
            for node in activenodes
                GLARM = globalarms(R, node) 

                preimages = filter!(x->x!==chpoint,findall(x->x==node,F)) #the characteristic point is the only one which has a danger of going twice
                append!(newactivenodes,preimages)

                for point in preimages
                    preimorder = Int[]
                    tuplelist = []


                    localarms = E[point]
                    for vertex in localarms
                        target = F[vertex]
                        x = findthe(target,GLARM)
                        push!(tuplelist,(x,vertex))
                    end

                    #sort according to first element of the tuple, the preimage order
                    sort!(tuplelist)
                    preimorder = [tuple[2] for tuple in tuplelist]

                    R[point] = preimorder

                end
            end
            activenodes = newactivenodes
        end
       

    end

    #We now want to add the critical orbit to R
    for ii in 1:z+4
        R[ii] = [x for x in E[ii]]
    end

    #now cyclicly permute this order so the arm towards zero is last
    for (point,arms) in enumerate(R)

        if point !==z && length(arms) > 1
            SECONDGLARM = globalarms(R,point)

            y = findthe(z,SECONDGLARM)
            circshift!(R[point],-y)
        end

    end

    return (R,F,markedpoints)
end

function findthe(item,listoflists)

    instances = [findall(x->x==item,list) for list in listoflists]

    hits = findall(x->!isempty(x),instances)

    if isempty(hits)
        error("no $item found")
    elseif length(hits) !== 1
        error("item $item found in multiple lists")
    end

    if length(instances[hits[1]]) == 1 
        return hits[1]        
    else
        v = hits[1]
        error("$item found more than once in list $v")
    end
end

function globalarms(E,point)
    neighbors = E[point]
    GA = [[point,n] for n in neighbors]
    for arm in GA
        activenodes = [arm[2]] 
        while !isempty(activenodes)
            newactivenodes = []
            for node in activenodes
                append!(newactivenodes,filter(x->!(x in arm),E[node]))
            end
            append!(arm,newactivenodes)
            activenodes = newactivenodes
        end
    end
    return GA
end

function branchpoints(E)
    Bindices = Int[]
    B = Vector{Int}[]
    for (ii,neighbors) in enumerate(E)
        if length(neighbors) > 2
            push!(B,neighbors)
            push!(Bindices)
        end
    end
    return B, Bindices
end

function arbembed((E,F,markedpoints))
    R = [collect(x) for x in E]
    z = findall(x->x==1,F)[1]
    for (point,arms) in enumerate(R)

        if point !==z && length(arms) > 1
            SECONDGLARM = globalarms(R,point)

            y = findthe(z,SECONDGLARM)
            circshift!(R[point],-y)
        end
    end
    return R
end


function characteristicpoints((E,F,markedpoints))

    P = Int[]

    for orb in components(F)
        c = characteristic(collect(orb),(E,F,markedpoints))
        if !isnothing(c)
            push!(P,c)
        end
    end

    return P
end

function characteristic(orb,(E,F,markedpoints)) #There is a unique characteristic point for each periodic orbit
    z = findall(x->x==1,F)[1]

    beta = z + 1
    mibeta = z + 2
    pa = z + 3
    a = z + 4
    if (1 in orb)
        return 
    elseif (beta in orb) 
        return 
    elseif (a in orb)
        return 
    else
        R = arbembed((E,F,markedpoints))
        #then there is a characteristic point in this orbit
        periodicpoints = filter(x->markedpoints[x].preperiod == 0,orb)
        for p in periodicpoints
            GLARM = globalarms(R,p)
            #All the other points on this orbit and zero are in one arm, and the critical point is in another
            #4.1 page 34 TreesBook
            targets = filter(x->x!==p,orb)
            push!(targets,z)

            arms = [findthe(point,GLARM) for point in targets]

            if allequal(arms)
                if findthe(1,GLARM) !== arms[1]
                    return p
                end
            end
        end
    end
end

function characteristicpoints(theta::Rational)
    return characteristicpoints(hubbardtree(theta))
end

function characteristicpoints(intadd::Vector{Int})
    return characteristicpoints(hubbardtree(intadd))
end

function binary(AIA::AngledInternalAddress)
    (R,F,markedpoints) = embedtree(AIA)
    digits = Char[]   

    z = findall(x->x==1,F)[1]
    pz = findall(x->x==z,F)[1]
    
    beta = z + 1
    mibeta = z + 2
    pa = z + 3
    a = z + 4

    dinger = 1

    for ii in 1:z
        x = markedpoints[ii].items[1]
        if x == '*'
            if dinger == 1
                push!(digits,'0')
            else
                push!(digits,'1')
            end
        elseif x == 'A'
            parent = R[ii][end]
            PARMS = globalarms(R,parent)
            b = findthe(mibeta,PARMS)
            point = findthe(ii,PARMS)
            if point == b
                if dinger == 1
                    push!(digits,'0')
                else
                    push!(digits,'1')
                end
            elseif point < b
                push!(digits,'0')
            else
                push!(digits,'1')
            end
        elseif x == 'B'
            parent = R[ii][end]
            PARMS = globalarms(R,parent)
            b = findthe(beta,PARMS)
            point = findthe(ii,PARMS)
            if point == b
                if dinger == 1
                    push!(digits,'1')
                else
                    push!(digits,'0')
                end
            elseif point < b
                push!(digits,'1')
            else
                push!(digits,'0')
            end
        else
            error("$x is in a kneading sequence!")
        end

        if markedpoints[pz].items[1] == 'A'
            (type,seq) = iteratetriod(markedpoints[1],(markedpoints[ii],markedpoints[pz],markedpoints[pa]))
        else
            (type,seq) = iteratetriod(markedpoints[1],(markedpoints[ii],markedpoints[pz],markedpoints[a]))
        end
        if type == "flat" && seq == 1  
            dinger = dinger*-1
            #println("$ii dinged")
        end

    end

    return Sequence(digits,0)
end

function test(intadd::Vector{Int})
    H = hubbardtree(intadd)
    n = length(characteristicpoints(H))
    num = fill(1,n)
    return internaladdress(angleof(binary(H,num)))
end

function adjlist(A::Matrix)
    n = size(A)[1]
    E = [Set{Int}() for kk in 1:n]
    for ii in 1:n
        for jj in 1:n
            if A[ii,jj] !== 0
                push!(E[ii],jj)
            end
        end
    end

    return E
end

function components(F::Vector{Int})
    C = Set{Int}[]
    done = zeros(Int64,length(F))

    while prod(done) == 0
        seedindex = findfirst(x->x==0,done)
        push!(C,Set(seedindex))
        activenodes = [seedindex]
        done[seedindex] = 1
        while !isempty(activenodes)
            newactivenodes = []
            for node in activenodes
                append!(newactivenodes,filter(x->!in(x,C[end]),append!(findall(x->x==node,F),F[node])))
            end
            for node in newactivenodes
                push!(C[end],node)
            end
            for index in findall(x->x in newactivenodes,1:length(F))
                done[index] = 1
            end
            activenodes = newactivenodes
        end
    end

    return C
end

