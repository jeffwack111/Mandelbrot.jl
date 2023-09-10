include("Sequences.jl")

function orbit(angle::Rational)

    events = Rational[]
    mapsto = Int[]

    while isempty(findall(x->x==angle,events))
        push!(events,angle)
        angle = angle*2
        angle = angle%1//1
    end

    for theta in events 
        push!(mapsto,findall(x->x==theta*2%1,events)[1])
    end

    return (events,mapsto)
    
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