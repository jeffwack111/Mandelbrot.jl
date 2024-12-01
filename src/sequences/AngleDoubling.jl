include("Sequences.jl")

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

abstract type KneadingSymbol end

struct A <: KneadingSymbol end
struct B <: KneadingSymbol end
struct star <: KneadingSymbol end

Base.show(io::IO, symb::A) = print(io,'A')
Base.show(io::IO, symb::B) = print(io,'B')
Base.show(io::IO, symb::star) = print(io,'*')

const KneadingSequence = Sequence{KneadingSymbol}

function KneadingSequence(angle::Rational)
    orb = orbit(angle)
    return thetaitinerary(angle,orb)
end

function thetaitinerary(theta::Rational,angle::Rational)
    return thetaitinerary(theta,orbit(angle))
end

function thetaitinerary(theta::Rational,orb::Sequence)
    a = theta/2
    b = (theta+1)/2
    itinerary = KneadingSymbol[]

    for angle in orb.items
        if angle == a
            push!(itinerary,star())
        elseif angle == b
            push!(itinerary,star())
        elseif angle > a && angle < b
            push!(itinerary,A())
        else
            push!(itinerary,B())
        end
    end
    
    return Sequence{KneadingSymbol}(itinerary,orb.preperiod)
end

function Base.show(io::IO, K::KneadingSequence)
    str = join([repr(item) for item in K.items])
    L = K.preperiod
    if L == 0
        return println(io,"|"*str*"|")
    else
        return println(io,str[1:L]*"|"*str[L+1:end]*"|")
    end
end

struct InternalAddress
    addr::Vector{Int}
end

function Base.length(intadd::InternalAddress)
    return length(intadd.addr)
end

function Base.getindex(intadd::InternalAddress,n::Int)
    return intadd.addr[n]
end

#TreesBook pg 21 def 2.2
function rho(K::KneadingSequence)
    function rho_nu(n::Int)
        if K.preperiod == 0 && mod(n,K.period) == 0
            return Inf
        end
        q = n+1 
        while K[q] == K[q-n]
            q = q+1
        end
        return q
    end
    return rho_nu
end

#works only for strictly periodic kneading sequences
function rhosequence(K::KneadingSequence,n::Int)
    if K.preperiod == 0
        seq = Vector{Int}([n])
        rho_nu = rho(K)
        while rho_nu(seq[end]) !== Inf
            append!(seq,rho_nu(seq[end]))
        end
        return seq
    else
        error("cant deal with preperiodic")
    end
end

function InternalAddress(K::KneadingSequence)
    return InternalAddress(rhosequence(K,1))
end

#Works only for finite internal addresses
function KneadingSequence(intadd::InternalAddress)
    address = copy(intadd.addr)
    if address == [1]
        return Sequence{KneadingSymbol}(KneadingSymbol[A()],0)
    else
        s = pop!(address)
        K = KneadingSequence(InternalAddress(address))
        R = K.items[mod1.(1:s-1,end)]
        if K.items[mod1(s,end)] == A()
            push!(R,B())
        else
            push!(R,A())
        end
        return Sequence{KneadingSymbol}(R,0)
    end
end

function admissible(K::KneadingSequence,m::Int)
    rho_nu = rho(K)
    if m in InternalAddress(K).addr
        return true
    else
        for k in 1:m-1
            if mod(m,k) == 0 && rho_nu(k) > m
                return true
            end
        end
        if rho_nu(m) != Inf
            for r in 1:m
                if mod(r,m) == mod(rho_nu(m),m) && ! (m in rhosequence(K,r))
                    return true
                end
            end
        end
        return false
    end
end

#Lemma 5.4 page 43 TreesBook
function admissible(K::KneadingSequence)
    for m in 1:K.period
        if !admissible(K,m)
            return (false,m)
        end
    end
    return true
end

function admissible(intadd::InternalAddress)
    return admissible(KneadingSequence(intadd))
end

struct AngledInternalAddress
    addr::Vector{Int}
    angles::Vector{Rational}
end

#Lemma 11.14 TreesBook page 146
function AngledInternalAddress(theta::Rational)
    intadd = InternalAddress(KneadingSequence(theta))
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

    return AngledInternalAddress(intadd.addr, angles)

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

#what is this for?
function orbit(a::Vector{Rational{Int}})
    items = Vector{Rational}[]
    while isempty(findall(x->x==a,items))
        push!(items,a)
        a = [(angle*2)%1//1 for angle in a]
    end

    preperiod = findall(x->x==a,items)[1] - 1

    return Sequence{Vector{Rational}}(items,preperiod)
    
end

function orbit(seq::Sequence)
    items = Sequence[]

    while isempty(findall(x->x==seq,items))
        push!(items,seq)
        seq = shift(seq)
    end

    preperiod = findall(x->x==seq,items)[1] - 1

    return Sequence{Sequence}(items,preperiod)
end

#TreesBook pg 145
function denominators(S::InternalAddress)
    n = length(S) - 1
    #the last entry has no angle #page 143

    K = KneadingSequence(S)

    q = Int[]

    for k in 1:n
        r = mod1(S[k+1],S[k])
        orb = rhosequence(K,r)
        if S[k] in orb
            qk = (S[k+1] - r)/S[k] + 1
        else
            qk = (S[k+1] - r)/S[k] + 2
        end
        push!(q,qk)
    end

    return q
end

function newdenominator(AIA::AngledInternalAddress,new::Int)
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
    return orbit(theta).period
end

function num_choices(D)
    return prod([d - 1 for d in D])
end

function num_angles(K)
    intadd = internaladdress(K)
    D = denominators(intadd)
    return num_choices(D)
end

function num_partners(theta)
    return num_angles(kneadingsequence(theta))
end

struct Digit{N} <: Integer
    value::Int
    function Digit{N}(value::Int) where N
        m = N-1
        0 <= value < N || error("Value must be in range [0, $m]")
        new{N}(value)
    end
end

Base.show(io::IO, d::Digit{N}) where {N} = print(io, d.value)
Base.Int(d::Digit{N}) where {N} = d.value

const BinaryExpansion = Sequence{Digit{2}}

function BinaryExpansion(theta::Rational)    
    orb = orbit(theta)
    itinerary = Digit{2}[]
    zero = Digit{2}(0)
    one = Digit{2}(1)
    for theta in orb.items
        if theta < 1//2
            push!(itinerary,zero)
        else
            push!(itinerary,one)
        end
    end
    return Sequence{Digit{2}}(collect(itinerary),orb.preperiod)
end

function Rational(bits::BinaryExpansion)
    theta = 0//1
    k = bits.period
    r = 1//(1//1-(2//1)^-k)
    for (ii,digit) in enumerate(bits.items)
        if digit==Digit{2}(1)
            if  ii<=bits.preperiod
                theta += (2//1)^(-ii)
            else 
                theta += r*(2//1)^(-ii)
            end
        end
    end
    return theta
end