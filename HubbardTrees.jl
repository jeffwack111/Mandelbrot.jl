using NetworkLayout
using IterTools

include("Sequences.jl")

#Definition 3.3
function triodmap(K::Sequence, arms::Tuple{Sequence,Sequence,Sequence})

    if arms[1].items[1] == arms[2].items[1] && arms[1].items[1] == arms[3].items[1]
        return (shift(arms[1]),shift(arms[2]),shift(arms[3]))

    elseif Set([arms[1].items[1],arms[2].items[1],arms[3].items[1]]) == Set(['A','B','*'])
        return "STOP"

    elseif arms[1].items[1] == arms[2].items[1] 
        return (shift(arms[1]),shift(arms[2]),K)

    elseif arms[1].items[1] == arms[3].items[1]
        return (shift(arms[1]),K,shift(arms[3]))

    elseif arms[2].items[1] == arms[3].items[1]
        return (K,shift(arms[2]),shift(arms[3]))

    else
        return error("no cases satisfied")

    end
end

function forwardorbit(S::Sequence)
    O = [S]
    for i in Iterators.drop(eachindex(S.items),1)
        push!(O,shift(O[end]))
    end
    return Sequence(O,S.pointsto)
end

function prependstar(S::Sequence)
    newpointsto = [2]
    newitems = ['*']
    append!(newitems,S.items)
    append!(newpointsto,S.pointsto .+ 1)
    return Sequence(newitems,newpointsto)
end

function S(K::Sequence)
    if K.pointsto[end] == 1
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
    for (s,t,u) in subsets(orb,3)

    end
end

function testbranch(K::Sequence,arms::Tuple{Sequence,Sequence,Sequence})
    triodsequence = Tuple{Sequence,Sequence,Sequence}[]
    triod = arms

    while isempty(findall(x->x==triod,triodsequence))
        push!(triodsequence,triod)
        triod = triodmap(K,triod)
    end

    mapsto = []

    for ii in Iterators.drop(eachindex(triodsequence),1)
        push!(mapsto,ii)
    end

    joint = findall(x->x==triod,triodsequence)[1]
    push!(mapsto,joint)

    return Sequence(triodsequence,mapsto)

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
    return Sequence(newitems,S.pointsto)
end
