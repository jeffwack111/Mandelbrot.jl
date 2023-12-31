using Primes
using IterTools

struct Sequence
    items::Vector
    preperiod::Int
end

import Base.==
function ==(x::Sequence,y::Sequence)
    return x.items==y.items && x.preperiod==y.preperiod
end

function Base.getindex(S::Sequence, ii::Int)
    k = length(S.items) - S.preperiod
    return S.items[mod1(ii-S.preperiod,k)]
end

function Base.getindex(S::Sequence, I::UnitRange)
    return [S[ii] for ii in I]
end

function period(S::Sequence)
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

function reduce(seq::Sequence)
    #first check that the periodic part is prime
    repetend = seq.items[(seq.preperiod+1):end]
    k = length(repetend)

    if length(repetend) > 1
        for d in divisors(k)
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
        return Sequence(append!(copy(preperiod),repetend),length(preperiod)) 
    else
        return Sequence(repetend,0)
    end

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

function divisors(n)
    d = Int64[1]
    for (p, e) in factor(n)
        r = 1
        l = length(d)
        for i in 1:e
            r *= p
            for j in 1:l
                push!(d, d[j]*r)
            end
        end
    end
    return sort(d)
end




