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

function binary(angle::Rational)
    
    orb = orbit(angle)
    itinerary = Char[]
    for theta in orb.items
        if theta < 1//2
            push!(itinerary,'0')
        else
            push!(itinerary,'1')
        end
    end
    return Sequence(collect(itinerary),orb.preperiod)

end

function internal_address(periodicsequence::Vector,mapsto::Vector{Int})
    if mapsto[end] == 1
        S = [1]
        v = ['A']
        for ii in eachindex(periodicsequence)
            if periodicsequence[ii] != v[mod1(ii,end)]
                push!(S,ii)
                v = periodicsequence[1:ii]
            end
        end
        return S
    else
        throw(DomainError("The internal address for a preperiodic sequence may be infinite"))
    end
end

function internal_address(angle::Rational)
    return internal_address(kneading_sequence(angle)...)    
end

function kneading_sequence(internaladdress::Vector{Int})
    if internal_address == [1]
        return PeriodicSequence([1])
    else
        s = pop!(internal_address)
        K = kneading_sequence(internal_address)
        R = K[1:s-1]
        if K[s] == 0
            push!(R,1)
        else
            push!(R,0)
        end
        return PeriodicSequence(R)
    end
end
