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

function goesto(S::Sequence)
    return push!(collect(2:length(S.items)),S.preperiod+1)
end

function shift(seq::Sequence)
    if seq.preperiod == 0 
        #then seq is periodic
        return Sequence(circshift!(copy(seq.items),-1),0)
    else
        #then the sequence is preperiodic
        return Sequence(collect(Iterators.drop(seq.items,1)),seq.preperiod-1)
    end  
end



    
    



