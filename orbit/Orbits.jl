using Primes
using IterTools

struct Orbit
    items::Vector
    preperiod::Int

    #A constructor must take a seed and a map. It better terminate! How do we make these safe to construct?

end

import Base.==
function ==(x::Orbit,y::Orbit)
    return x.items==y.items && x.preperiod==y.preperiod
end

function Base.getindex(S::Orbit, ii::Int)
    k = length(S.items) - S.preperiod
    return S.items[mod1(ii-S.preperiod,k)]
end

function Base.getindex(S::Orbit, I::UnitRange)
    return [S[ii] for ii in I]
end

function period(S::Orbit)
    return length(S.items) - S.preperiod
end

function removerepetend(r::Vector,p::Vector)
    k = length(r)
    if length(r) > length(p) #then p certainly does not contain r
        return p
    else
        if p[(end-k+1):end] == r
            removerepetend(r,p[1:(end-k)])
        else
            return p
        end
    end
end

function reduce(seq::Orbit)
    #first check that the periodic part is prime
    repetend = seq.items[(seq.preperiod+1):end]
    k = length(repetend)

    if length(repetend) > 1
        for d in Iterators.flatten((1,factor(Set,k)))
            chunks = collect.(partition(repetend,d))
            if allequal(chunks)
                repetend = chunks[1]
                k = length(repetend)
                break
            end
        end
    end

    if seq.preperiod != 0
        preperiod = copy(seq.items[1:seq.preperiod])
        while !isempty(preperiod) && repetend[end] == preperiod[end]
            pop!(preperiod)
            circshift!(repetend,1)
        end
        return Orbit(append!(copy(preperiod),repetend),length(preperiod)) 
    else
        return Orbit(repetend,0)
    end

end

    

function goesto(S::Orbit)
    return push!(collect(2:length(S.items)),S.preperiod+1)
end

function shift(seq::Orbit)
    if seq.preperiod == 0 
        #then seq is periodic
        return Orbit(circshift(copy(seq.items),-1),0)
    else
        #then the sequence is preperiodic
        return Orbit(collect(Iterators.drop(seq.items,1)),seq.preperiod-1)
    end  
end



    
    



