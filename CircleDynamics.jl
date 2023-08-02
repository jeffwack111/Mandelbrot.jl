function orbit(angle::Rational)

    orb = Rational[]

    while isempty(findall(x->x==angle,orb))
        push!(orb,angle)
        angle = angle*2
        angle = angle%1//1
    end

    l = findall(x->x==angle,orb)[1]
    if l == 1
        return PVector(orb)
    else 
        return PPVector(orb[1:(l-1)],orb[l:end])
    end

end

function binary(angle::Rational)

    if angle >= 1 || angle < 0
        throw(DomainError(angle,"an angle is measured from 0 to 1 full turn"))
    end

    self_itinerary = Char[]
    orb = []

    while isempty(findall(x->x==angle,orb))
        push!(orb,angle)
        if angle < 1//2
            push!(self_itinerary,'0')
        else
            push!(self_itinerary,'1')
        end
        angle = angle*2
        angle = angle%1//1
    end

    l = findall(x->x==angle,orb)[1]
    if l == 1
        return PString(join(self_itinerary))
    else 
        return PPString(join(self_itinerary[1:(l-1)]),join(self_itinerary[l:end]))
    end

end

function kneading_sequence(angle::Rational)
    
    a = angle/2
    b = (angle+1)/2
    self_itinerary = Char[]
    orb = []

    while isempty(findall(x->x==angle,orb))
        push!(orb,angle)
        if angle == a
            push!(self_itinerary,'a')
        elseif angle == b
            push!(self_itinerary,'b')
        elseif angle > a && angle < b
            push!(self_itinerary,'A')
        else
            push!(self_itinerary,'B')
        end
        angle = angle*2
        angle = angle%1//1
    end

    l = findall(x->x==angle,orb)[1]
    if l == 1
        return PString(join(self_itinerary))
    else 
        return PPString(join(self_itinerary[1:(l-1)]),join(self_itinerary[l:end]))
    end

end