include("Sequences.jl")

struct AngledInternalAddress
    addr::Vector{Int}
    angles::Vector{Rational}
end

#Lemma 11.14 TreesBook page 146
function AngledInternalAddress(theta::Rational)
    intadd = internaladdress(theta)
    denoms = denominators(intadd)
    angles = Rational[]

    for ii in 1:length(intadd)-1
        n = 0
        for jj in 1:(denoms[ii] - 1)
            if ((2//1)^((jj-1)*intadd[ii])*theta)%1 <= theta
                n += 1
            end
        end
        push!(angles,n//denoms[ii])
    end

    return AngledInternalAddress(intadd, angles)

end

function firstaddress(intadd::Vector{Int})
    d = denominators(intadd)
    angles = [1//x for x in d]
    return AngledInternalAddress(intadd,angles)
end

function bifurcate(aia::AngledInternalAddress)
    intadd = copy(aia.addr)
    intadd = append!(intadd,intadd[end]*2//1)
    angles = append!(copy(aia.angles),1//2)
    return AngledInternalAddress(intadd, angles)
end


function nextaddress(AIA::AngledInternalAddress)
    #give the next valid internal address, going counterclockwise
end

function orbit(angle::Rational)
    items = Rational[]

    while isempty(findall(x->x==angle,items))
        push!(items,angle)
        angle = angle*2
        angle = angle%1//1
    end

    preperiod = findall(x->x==angle,items)[1] - 1

    return Sequence{Rational}(items,preperiod)
    
end

function orbit(a::Vector{Rational{Int}})
    items = Vector{Rational}[]
    while isempty(findall(x->x==a,items))
        push!(items,a)
        a = [(angle*2)%1//1 for angle in a]
    end

    preperiod = findall(x->x==a,items)[1] - 1

    return Sequence{Vector{Rational}}(items,preperiod)
    
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
    
    return Sequence{Char}(itinerary,orb.preperiod)
end

function thetaitinerary(theta::Rational,angle::Rational)
    return thetaitinerary(theta,orbit(angle))
end

function kneadingsequence(angle::Rational)
    orb = orbit(angle)
    return thetaitinerary(angle,orb)
end

function binary(theta::Rational)
    
    orb = orbit(theta)
    itinerary = Char[]
    for theta in orb.items
        if theta < 1//2
            push!(itinerary,'0')
        else
            push!(itinerary,'1')
        end
    end
    return Sequence{Char}(collect(itinerary),orb.preperiod)

end


#this better be a sequence of 1s and 0s
function angleof(binary::Sequence)
    theta = 0//1
    k = period(binary)
    r = 1//(1//1-(2//1)^-k)
    for (ii,digit) in enumerate(binary.items)
        if digit=='1' 
            if  ii<=binary.preperiod
                theta += (2//1)^(-ii)
            else 
                theta += r*(2//1)^(-ii)
            end
        end
    end
    return theta
end

function angleof(digs::String)
    return angleof(Sequence{Char}(collect(digs),0))
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

function rhoSequence(kneadingsequence::Sequence,n::Int) 
    if kneadingsequence.preperiod == 0
        seq = Int[n]
        while rho(kneadingsequence,seq[end]) !== Inf
            append!(seq,rho(kneadingsequence,seq[end]))
        end
        return seq
    else
        error("cant deal with preperiodic")
    end
end

function rhoSequence(K::Sequence)
    function r(n::Int)
        return rhoSequence(K,n)
    end
    return r
end

#works only for strictly periodic kneading sequences
function internaladdress(K::Sequence)
    return rhoSequence(K,1)
end

function internaladdress(theta::Rational)
    return rhoSequence(kneadingsequence(theta),1)    
end

#Works only for finite internal addresses
function kneadingsequence(intadd::Vector{Int})
    internaladdress = copy(intadd)
    if internaladdress == [1]
        return Sequence{Char}(['A'],0)
    else
        s = pop!(internaladdress)
        K = kneadingsequence(internaladdress)
        R = K.items[mod1.(1:s-1,end)]
        if K.items[mod1(s,end)] == 'A'
            push!(R,'B')
        else
            push!(R,'A')
        end
        return Sequence{Char}(R,0)
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
                if mod(r,m) == mod(rho(kneadingsequence, m),m) && ! (m in rhoSequence(kneadingsequence,r))
                    return true
                end
            end
        end
        return false
    end
end

#Lemma 5.4 page 43 TreesBook
function admissible(K::Sequence)
    for m in 1:period(K)
        if !admissible(K,m)
            return (false,m)
        end
    end
    return true
end

function admissible(intadd::Vector{Int})
    return admissible(kneadingsequence(intadd))
end

#TreesBook pg 145
function denominators(S::Vector{Int})
    n = length(S) - 1
    #the last entry has no angle #page 143

    K = kneadingsequence(S)
    p = rhoSequence(K)

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

function denominator(AIA::AngledInternalAddress,new::Int)
    K = kneadingsequence(AIA.addr)
    p = rhoSequence(K)

    r = mod1(new,AIA.addr[end])
    orb = p(r)

    if AIA.addr[end] in orb
        qk = Int((new - r)/AIA.addr[end] + 1) #why is this always an integer?

    else
        qk = Int((new - r)/AIA.addr[end] + 2)

    end
    
    return qk
end


function numembeddings(I::Vector{Int})
    D = denominators(I)
    n = 1
    for d in D
        n = n*totient(d) #The number of integers less than d and coprime
    end
    return n
end

function numembeddings(theta::Rational)
    return numembeddings(internaladdress(theta))
end

function period(theta::Rational)
    return period(orbit(theta))
end



