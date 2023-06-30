struct PeriodicSequence{T} <: AbstractVector{T}
    S::Vector{T}
end

function Base.getindex(P::PeriodicSequence,i::Int)
    return P.S[mod1(i,end)]
end

function Base.getindex(P::PeriodicSequence, I) 
    return [P[i] for i in I]
end

function Base.size(P::PeriodicSequence)
    return Base.size(P.S)
end

struct PrePeriodicSequence{T} <: AbstractVector{T}
    L::Vector{T}
    K::Vector{T}
end

function Base.size(P::PrePeriodicSequence)
    return Base.size(P.L).+Base.size(P.K)
end

function Base.getindex(P::PrePeriodicSequence,i::Int)
    l = length(P.L)
    if i<=0
        throw(DomainError(i, "index must be positive"))
    elseif i <= l
        return P.L[i]
    else
        return P.K[mod1(i - l,end)]
    end
end

function Base.getindex(P::PrePeriodicSequence, I) 
    return [P[i] for i in I]
end