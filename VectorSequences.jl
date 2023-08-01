struct PVector{T} 
    S::Vector{T}
end

function Base.getindex(P::PVector,i::Int)
    return P.S[mod1(i,end)]
end

function Base.getindex(P::PVector, I) 
    return [P[i] for i in I]
end

struct PPVector{T}
    L::Vector{T}
    K::PVector{T}
end

function PPVector(L::Vector,K::Vector)
    return PPVector(L,PVector(K))
end

function Base.getindex(P::PPVector,i::Int)
    if i<1
        throw(BoundsError(P,i))
    end
    try
        return P.L[i]
    catch
        return P.K[i]
    end
end

function Base.getindex(P::PPVector, I) 
    return [P[i] for i in I]
end