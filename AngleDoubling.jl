include("Sequences.jl")

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

#works only for strictly periodic kneading sequences
function internaladdress(K::Sequence)
    if K.preperiod == 0
        S = [1]
        v = ['A']
        for ii in eachindex(K.items)
            if K.items[ii] != v[mod1(ii,end)]
                push!(S,ii)
                v = K.items[mod1.(1:ii,end)]
            end
        end
        return S
    else
        throw(DomainError("The internal address for a preperiodic sequence may be infinite"))
    end
end

function internaladdress(angle::Rational)
    return internaladdress(kneadingsequence(angle))    
end

#Works only for finite internal addresses
function kneadingsequence(internaladdress::Vector{Int})
    if internaladdress == [1]
        return Sequence(['A'],0)
    else
        s = pop!(internaladdress)
        K = kneadingsequence(internaladdress)
        R = K.items[mod1.(1:s-1,end)]
        if K.items[mod1(s,end)] == 'A'
            push!(R,'B')
        else
            push!(R,'A')
        end
        return Sequence(R,0)
    end
end
