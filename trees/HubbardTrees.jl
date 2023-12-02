#Code written by Jeff Wack
#Algorithms from
#Admissibility of kneading sequences and structure of Hubbard trees for quadratic polynomials
#Henk Bruin, Dierk Schleicher

include("AngleDoubling.jl")

function hubbardtree(seq::Sequence)
    orbit = criticalorbit(seq)
    n = length(orbit)

    results = fill(0,(n,n,n))

    markedpoints = copy(orbit)

    for triple in subsets(collect(enumerate(orbit)),3)
        type,seq = iteratetriod(orbit[2],(triple[1][2],triple[2][2],triple[3][2]))

        if type == "branched"
            if isempty(findall(x->x==seq,markedpoints))
                push!(markedpoints,seq)
                results = cat(results,1,dims=(1,2,3))
                index = lastindex(markedpoints)
            else
                index = findall(x->x==seq,markedpoints)[1]
            end
        elseif type == "flat"
            index = triple[seq][1]
        end

        results[triple[1][1],triple[2][1],triple[3][1]] = index

    end

    for triple in subsets(collect(enumerate(markedpoints)),3)
        if results[triple[1][1],triple[2][1],triple[3][1]] == 0 
            type,seq = iteratetriod(orbit[2],(triple[1][2],triple[2][2],triple[3][2]))

            if type == "branched"
                index = findall(x->x==seq,markedpoints)[1]
            elseif type == "flat"
                index = triple[seq][1]
            end

            results[triple[1][1],triple[2][1],triple[3][1]] = index
        end
    end

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

    A = zeros(Int64,(d,d))

    for (ii,jj) in subsets(collect(1:d),2)
        pointsbetween = 0
        for mid in M[ii,jj]
            if mid != ii && mid !=jj
                pointsbetween += 1
            end
        end
        if pointsbetween == 0
            A[ii,jj] = 1
            A[jj,ii] = 1
        end
    end

    return A, markedpoints
                    

end

function hubbardtree(angle::Rational)
    return hubbardtree(kneadingsequence(angle))
end

function hubbardtree(internaladdress::Vector{Int})
    K = kneadingsequence(internaladdress)
    seq = copy(K.items)
    seq[end] = '*'
    print(seq)
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
    return reduce(Sequence(newitems,S.preperiod))
end

function prependstar(S::Sequence)
    if S.preperiod == 0
        return Sequence(circshift!(copy(S.items),1),0)
    else
        return Sequence(pushfirst!(copy(S.items),'*'),S.preperiod+1)
    end
end

#Def. page 259
function criticalorbit(K::Sequence)
    if K.preperiod == 0
        A = [K]
        for i in Iterators.drop(eachindex(K.items),1)
            push!(A,shift(A[end]))
        end
        circshift!(A,1)
        return A  
    else
        A = [prependstar(K)]
        for i in eachindex(K.items)
            push!(A,shift(A[end]))
        end
        return A    
    end
end

function S(angle::Rational)
    return S(kneadingsequence(angle))
end

function V(K::Sequence)
    orb = S(K)
    branchpoints = Sequence[]
    for (s,t,u) in subsets(orb,3)
        type,seq = iteratetriod(K,(s,t,u))
        if type == "branched" && isempty(findall(x->x==seq,branchpoints))
            push!(branchpoints,seq)
        end
    end
    return append!(orb,branchpoints)
end

function V(angle::Rational)
    return V(kneadingsequence(angle))
end

