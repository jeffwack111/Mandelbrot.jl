using Primes
using IterTools

struct Digit{N} <: Integer
    value::Int
    function Digit{N}(value::Int) where N
        m = N-1
        0 <= value < N || error("Value must be in range [0, $m]")
        new{N}(value)
    end
end

Base.show(io::IO, d::Digit{N}) where {N} = print(io, d.value)
Base.Int(d::Digit{N}) where {N} = d.value

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

const BinaryExpansion = Sequence{Digit{2}}

function BinaryExpansion(theta::Rational)    
    orb = orbit(theta)
    itinerary = Digit{2}[]
    zero = Digit{2}(0)
    one = Digit{2}(1)
    for theta in orb.items
        if theta < 1//2
            push!(itinerary,zero)
        else
            push!(itinerary,one)
        end
    end
    return Sequence{Digit{2}}(collect(itinerary),orb.preperiod)
end

function orbit(angle::Rational)
    items = Rational[]

    while isempty(findall(x->x==angle,items))
        push!(items,angle)
        angle = angle*2
        angle = angle%1//1
    end

    preperiod = findall(x->x==angle,items)[1] - 1

    return Sequence{Rational}(items,preperiod)
    
end

abstract type KneadingSymbol end

struct A <: KneadingSymbol end
struct B <: KneadingSymbol end
struct star <: KneadingSymbol end

Base.show(io::IO, symb::KneadingSymbol) = print(io, repr(typeof(symb)))

const KneadingSequence = Sequence{KneadingSymbol}

function thetaitinerary(theta::Rational,orb::Sequence)
    a = theta/2
    b = (theta+1)/2
    itinerary = KneadingSymbol[]

    for angle in orb.items
        if angle == a
            push!(itinerary,star())
        elseif angle == b
            push!(itinerary,star())
        elseif angle > a && angle < b
            push!(itinerary,A())
        else
            push!(itinerary,B())
        end
    end
    
    return Sequence{KneadingSymbol}(itinerary,orb.preperiod)
end

function thetaitinerary(theta::Rational,angle::Rational)
    return thetaitinerary(theta,orbit(angle))
end

function KneadingSequence(angle::Rational)
    orb = orbit(angle)
    return thetaitinerary(angle,orb)
end



function Sequence(str::String,preperiod::Int)
    return Sequence{Char}(Vector{Char}(str),preperiod)
end

import Base.==
function ==(x::Sequence,y::Sequence)
    return x.items==y.items && x.preperiod==y.preperiod
end

function Base.getindex(S::Sequence, ii::Int)
    if ii <= S.preperiod
        return S.items[ii]
    else
        k = length(S.items) - S.preperiod
        return S.items[mod1(ii-S.preperiod,k)]
    end
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

function shiftby(seq::Sequence,n::Int)
    if n<0
        error("cannot shift a sequence a negative number of times")
    elseif n==0
        return seq
    elseif n==1
        return shift(seq)
    else
        return shift(shiftby(seq,n-1))
    end
end

function prepend(K::Sequence,thing)
    return Sequence{Char}(pushfirst!(copy(K.items),thing),K.preperiod+1)
end

