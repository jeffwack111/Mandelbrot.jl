#Admissibility of kneading sequences and structure of Hubbard trees for quadratic polynomials
#Henk Bruin, Dierk Schleicher

using NetworkLayout
using IterTools

include("AngleDoubling.jl")

#Definition 3.3, page 260
function triodmap(K::Sequence, arms::Tuple{Sequence,Sequence,Sequence})

    if arms[1].items[1] == arms[2].items[1] && arms[1].items[1] == arms[3].items[1]
        return (shift(arms[1]),shift(arms[2]),shift(arms[3])),0

    elseif Set([arms[1].items[1],arms[2].items[1],arms[3].items[1]]) == Set(['A','B','*'])
        return "STOP",0

    elseif arms[1].items[1] == arms[2].items[1] 
        return (shift(arms[1]),shift(arms[2]),K),3

    elseif arms[1].items[1] == arms[3].items[1]
        return (shift(arms[1]),K,shift(arms[3])),2

    elseif arms[2].items[1] == arms[3].items[1]
        return (K,shift(arms[2]),shift(arms[3])),1

    else
        return error("no cases satisfied")

    end
end

function forwardorbit(S::Sequence)
    O = [S]
    for i in Iterators.drop(eachindex(S.items),1)
        push!(O,shift(O[end]))
    end
    return Sequence(O,S.preperiod)
end

function prependstar(S::Sequence)
    if S.preperiod == 0
        return Sequence(circshift!(copy(S.items),1),0)
    else
        return Sequence(pushfirst!(copy(S.items),'*'),S.preperiod+1)
    end
end

#Def. page 259
function S(K::Sequence)
    if K.preperiod == 0
        A = [K]
        for i in Iterators.drop(eachindex(K.items),1)
            push!(A,shift(A[end]))
        end
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
        type,seq = testbranch(K,(s,t,u))
        if type == "branched" && isempty(findall(x->x==seq,branchpoints))
            push!(branchpoints,seq)
        end
    end
    return append!(orb,branchpoints)
end

function V(angle::Rational)
    return V(kneadingsequence(angle))
end

function testbranch(K::Sequence,arms::Tuple{Sequence,Sequence,Sequence})
    chopsequence = []
    chop = 0
    triodsequence = Tuple{Sequence,Sequence,Sequence}[]
    triod = arms

    while isempty(findall(x->x==triod,triodsequence)) && triod != "STOP"
        push!(triodsequence,triod)
        push!(chopsequence,chop) 
        triod,chop = triodmap(K,triod)
    end

    if triod == "STOP"
        return "flat"
    end

    preperiod = findall(x->x==triod,triodsequence)[1] - 1

    if (1 in chopsequence) && (2 in chopsequence) && (3 in chopsequence) 
        return "branched",majorityvote(Sequence(triodsequence,preperiod))
    else
        return "flat"
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
