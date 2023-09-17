struct Sequence
    items::Vector
    preperiod::Int
end

import Base.==
function ==(x::Sequence,y::Sequence)
    return x.items==y.items && x.preperiod==x.preperiod
end

function goesto(S::Sequence)
    return push!(collect(2:length(S.items)),S.preperiod+1)
end


function shift(seq::Sequence)
    if seq.preperiod == 0 
        #then seq is periodic

        newitems = circshift(items,-1)
        return Sequence(newitems,0)

    else
        #then the sequence is preperiodic

        newpreperiod = preperiod - 1 
        newitems = collect(Iterators.drop(seq.items,1))
        return Sequence(newitems,newpreperiod)

    end  
end

function orbit(angle::Rational)
    items = Rational[]

    while isempty(findall(x->x==angle,items))
        push!(items,angle)
        angle = angle*2
        angle = angle%1//1
    end

    preperiod = findall(x->x==angle,items)[1] - 1

    return Sequence(items,preperiod)
    
end

function kneadingsequence(orb::Sequence)
    angle = orb.items[1]

    a = angle/2
    b = (angle+1)/2
    self_itinerary = Char[]

    for theta in orb.items
        if theta == a
            push!(self_itinerary,'*')
        elseif theta == b
            push!(self_itinerary,'*')
        elseif theta > a && theta < b
            push!(self_itinerary,'A')
        else
            push!(self_itinerary,'B')
        end
    end
    
    return Sequence(self_itinerary,orb.preperiod)
end

function kneadingsequence(angle::Rational)
    orb = orbit(angle)
    return kneadingsequence(orb)
end
    
    



