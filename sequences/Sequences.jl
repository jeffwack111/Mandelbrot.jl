using Primes
using IterTools

struct Sequence{T}
    items::Vector{T}
    preperiod::Int

    #this inner constructor ensures that all sequences are in reduced form
    function Sequence{T}(items::Vector{T},preperiod::Int) where T
        #first check that the periodic part is prime
        repetend = items[(preperiod+1):end]
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
    
        if preperiod != 0
            preperiod = copy(items[1:preperiod])
            while !isempty(preperiod) && repetend[end] == preperiod[end]
                pop!(preperiod)
                circshift!(repetend,1)
            end
            return new{T}(append!(copy(preperiod),repetend),length(preperiod)) 
        else
            return new{T}(repetend,0)
        end
    
    end
end

function Sequence(str::String,preperiod::Int)
    return Sequence(Vector{Char}(str),preperiod)
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

function Base.length(S::Sequence)
    return length(S.items)
end

function Base.hash(S::Sequence,h::UInt)
    k = hash(S.preperiod,foldr(hash,S.items; init = UInt(0)))
    return hash(h,hash(:Sequence,k))
end

function Base.show(io::IO, s::Sequence{Char})
    str = String(s.items)
    L = s.preperiod
    if L == 0
        return print(io,"|"*str*"|")
    else
        return print(io,str[1:L]*"|"*str[L+1:end]*"|")
    end
end

function period(S::Sequence)
    return length(S.items) - S.preperiod
end

function shift(seq::Sequence{T}) where T
    if seq.preperiod == 0 
        #then seq is periodic
        return Sequence{T}(circshift(copy(seq.items),-1),0)
    else
        #then the sequence is preperiodic
        return Sequence{T}(collect(Iterators.drop(seq.items,1)),seq.preperiod-1)
    end  
end

function orbit(seq::Sequence)
    items = Sequence[]

    while isempty(findall(x->x==seq,items))
        push!(items,seq)
        seq = shift(seq)
    end

    preperiod = findall(x->x==seq,items)[1] - 1

    return Sequence{Sequence}(items,preperiod)
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




