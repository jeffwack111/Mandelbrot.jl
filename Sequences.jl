struct PeriodicSequence{T} <: AbstractVector{T}
    S::Vector{T}
end

function Base.getindex(P::PeriodicSequence,i::Int)
    return P.S[mod1(i,end)]
end

function Base.size(P::PeriodicSequence)
    return Base.size(P.S)
end

function Base.getindex(P::PeriodicSequence, I) 
    return [P[i] for i in I]
end