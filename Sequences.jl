#############################################

struct PString
    S::String
end

function Base.getindex(S::PString,i::Int)
    return S.S[mod1(i,end)]
end

function Base.getindex(P::PString, I) 
    return join([P[i] for i in I])
end

struct PPString
    L::String
    K::PString
end

function Base.getindex(S::PPString,i::Int)
    if i<1
        throw(BoundsError(S,i))
    end
    try
        return S.L[i]
    catch
        return S.K[i]
    end
end

function Base.getindex(S::PPString, I) 
    return join([S[i] for i in I])
end

function PPString(L::String, K::String)
    return PPString(L,PString(K))
end



############################################

struct PVector{T} 
    S::Vector{T}
end

function Base.getindex(P::PVector,i::Int)
    return P.S[mod1(i,end)]
end

function Base.getindex(P::PVector, I) 
    return [P[i] for i in I]
end

function Base.length(P::PVector)
    return length(P.S)
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

function Base.length(P::PPVector)
    return length(P.L) + length(P.K)
end