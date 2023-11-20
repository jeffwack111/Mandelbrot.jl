include("Sequences.jl")

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

function thetaitinerary(theta::Rational,orb::Sequence)
    a = theta/2
    b = (theta+1)/2
    itinerary = Char[]

    for angle in orb.items
        if angle == a
            push!(itinerary,'*')
        elseif angle == b
            push!(itinerary,'*')
        elseif angle > a && angle < b
            push!(itinerary,'A')
        else
            push!(itinerary,'B')
        end
    end
    
    return Sequence(itinerary,orb.preperiod)
end

function kneadingsequence(angle::Rational)
    orb = orbit(angle)
    return thetaitinerary(angle,orb)
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


function rho(kneadingsequence::Sequence,n::Int)
    if kneadingsequence.preperiod == 0 && mod(n,period(kneadingsequence)) == 0
        return Inf
    end
    k = n+1 
    while kneadingsequence[k] == kneadingsequence[k-n]
        k = k+1
    end
    return k
end

function rhoorbit(kneadingsequence::Sequence,n::Int) 
    if kneadingsequence.preperiod == 0
        orbit = [n]
        while rho(kneadingsequence,orbit[end]) !== Inf
            append!(orbit,rho(kneadingsequence,orbit[end]))
        end
        return orbit
    else
        error("cant deal with preperiodic")
    end
end

function rhoorbit(K::Sequence)
    function r(n::Int)
        return rhoorbit(K,n)
    end
    return r
end

#works only for strictly periodic kneading sequences
function internaladdress(K::Sequence)
    return rhoorbit(K,1)
end

function internaladdress(angle::Rational)
    return rhoorbit(kneadingsequence(angle),1)    
end

#Works only for finite internal addresses
function kneadingsequence(intadd::Vector{Int})
    internaladdress = copy(intadd)
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

function admissible(kneadingsequence::Sequence,m::Int)
    if m in internaladdress(kneadingsequence)
        return true
    else
        for k in 1:m-1
            if mod(m,k) == 0 && rho(kneadingsequence,k) > m
                return true
            end
        end
        if rho(kneadingsequence, m) != Inf
            for r in 1:m
                if mod(r,m) == mod(rho(kneadingsequence, m),m) && ! (m in rhoorbit(kneadingsequence,r))
                    return true
                end
            end
        end
        return false
    end
end

#TreesBook pg 145
function denominators(S::Vector{Int})
    n = length(S) - 1
    #the last entry has no angle #page 143

    K = kneadingsequence(S)
    p = rhoorbit(K)

    q = Int[]

    for k in 1:n
        r = mod1(S[k+1],S[k])
        orb = p(r)
        if S[k] in orb
            qk = (S[k+1] - r)/S[k] + 1

        else
            qk = (S[k+1] - r)/S[k] + 2

        end
        push!(q,qk)
    end

    return q
end


