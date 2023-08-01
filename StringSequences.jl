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