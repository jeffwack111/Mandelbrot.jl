using Primes

struct Sequence
    items::Vector
    preperiod::Int
end

import Base.==
function ==(x::Sequence,y::Sequence)
    return x.items==y.items && x.preperiod==x.preperiod
end

function Base.getindex(S::Sequence, ii::Int)
    k = length(S.items) - S.preperiod
    return S.items[mod1(ii-S.preperiod,k)]
end

function Base.getindex(S::Sequence, I::UnitRange)
    return [S[ii] for ii in I]
end

function removerepetend(r::Vector,p::Vector)
    k = length(r)
    if length(r) > length(p) #then p certainly does not contain ratio
        return p
    else
        if p[(end-k+1):end] == r
            removerepetend(r,p[1:(end-k)])
        else
            return p
        end
    end
end

function reduce(seq::Sequence)
    #first check that the periodic part is prime
    if seq.preperiod != 0
        preperiod = seq.items[1:seq.preperiod]
    else
        preperiod = empty(seq.items)
    end

    repetend = seq.items[(seq.preperiod+1):end]
    k = length(repetend)

    for d in Iterators.flatten((1,factor(Set,k)))
        chunks = collect.(partition(repetend,d))
        if allequal(chunks)
            repetend = chunks[1]
            k = length(repetend)
            break
        end
    end

    preperiod = removerepetend(repetend,preperiod)
    l = length(preperiod)

    return Sequence(append!(preperiod,repetend),l) 
end

    

function goesto(S::Sequence)
    return push!(collect(2:length(S.items)),S.preperiod+1)
end

function shift(seq::Sequence)
    if seq.preperiod == 0 
        #then seq is periodic
        return Sequence(circshift(copy(seq.items),-1),0)
    else
        #then the sequence is preperiodic
        return Sequence(collect(Iterators.drop(seq.items,1)),seq.preperiod-1)
    end  
end



    
    



