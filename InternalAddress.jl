function kneading_sequence(angle::Rational)

    a = angle/2
    b = (angle+1)/2
    self_itinerary = Char[]
    orb,mapsto = orbit(angle)

    for theta in orb
        if theta == a
            push!(self_itinerary,'a')
        elseif theta == b
            push!(self_itinerary,'b')
        elseif theta > a && theta < b
            push!(self_itinerary,'A')
        else
            push!(self_itinerary,'B')
        end
    end

    return self_itinerary,mapsto
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
    if A == [1]
        return PeriodicSequence([1])
    else
        s = pop!(A)
        K = kneading_sequence(A)
        R = K[1:s-1]
        if K[s] == 0
            push!(R,1)
        else
            push!(R,0)
        end
        return PeriodicSequence(R)
    end
end
