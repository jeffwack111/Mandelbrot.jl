using IterTools
include("SpiderFuncs.jl")
include("../sequences/AngleDoubling.jl")

struct Spider
    angle::Rational
    orbit::Sequence{Rational}
    kneading_sequence::Sequence{Char}
    legs::Vector{Vector{ComplexF64}}
end

function Base.show(io::IO, S::Spider)
    npoints = sum([length(leg) for leg in S.legs])
    return println("A "*repr(S.angle)*"-spider with kneading sequence "*repr(S.kneading_sequence)*" and "*repr(npoints)*" points.")
end

function standardspider(theta::Rational)
    orb = orbit(theta)
    K = thetaitinerary(theta,orb) #the kneading sequence is calculated like so to avoid recalculating the orbit of theta
    legs = standardlegs(orb)
    return Spider(theta,orb,K,legs)
end

function standardlegs(orbit::Sequence)
    r = collect(LinRange(100,1,10))  #NOTE - it may be better to keep this as a 'linrange' but I don't understand what that means
    legs  = Vector{ComplexF64}[]
    for theta in orbit.items
        push!(legs,(cos(theta*2*pi)+1.0im*sin(theta*2*pi)) .* r)
    end
    legs[1] = legs[1] .- legs[1][end]
    return legs
end

function standardlegs(angle::Rational)
    info = SpiderInfo(angle)
    return standardlegs(info)
end

function grow!(legs::Vector{Vector{ComplexF64}},scale::Real,num::Int)
    #Get the arrow pointing from the second to last point to the last point
    for leg in legs
        arrow = scale*(leg[1]-leg[2])
        #add new elements which are the last element plus the arrow
        for i in 1:num
            pushfirst!(leg,leg[1]+arrow)
        end
    end
end


function prune!()
end

function stats(legs::Vector{Vector{ComplexF64}})
    println("~~~~~~~~~~~~~~~~~~~~")
    λ = legs[2][end]
    println(λ/2)
    points = length(leg)*length(leg[1])
    println(points)
    radius = sqrt(abs2(leg[1][1]))
    println(radius)
end

function mapspider!(S::Spider)
    λ = S.legs[2][end] #the parameter of our polynomial map z -> λ(1+z/2)^2
    a = (angle(S.legs[1][1]/λ)/(2*pi)+1)%1 #the angle whose halves will define the boundary between regions A and B at infinity. in full turns

    n = length(S.orbit)

    newLegs = Vector{Vector{ComplexF64}}(undef,n)

    #leg 1 goes first
    newLegs[1] = path_sqrt(S.legs[2]./λ)

    theta1 = (angle(newLegs[1][1])/(2*pi)+1)%1 #gives the angle of the shoulder in full turns
    #this angle lays in region A by definition. 

    if a/2 < theta1 && theta1 < (a+1)/2
        cregion = 'A'
        dregion = 'B'
    else
        cregion = 'B'
        dregion = 'A'
    end
    
    for ii in 2:n-1
        #first find the right half plane preimage of the shoulder
        u = sqrt(S.legs[ii+1][1]./λ)
        thetau = (angle(u)/(2*pi)+1)%1

        if a/2 < thetau && thetau < (a+1)/2 #thetau is in the connected region
            if S.kneading_sequence[ii] == cregion #then this is the correct preimage
                newLegs[ii] = path_sqrt(S.legs[ii+1]./λ)
            else
                newLegs[ii] = -1 .* path_sqrt(S.legs[ii+1]./λ)
            end
        else #thetau is in the disconnected region
            if S.kneading_sequence[ii] == dregion #then this is the correct preimage
                newLegs[ii] = path_sqrt(S.legs[ii+1]./λ)
            else
                newLegs[ii] =  -1 .* path_sqrt(S.legs[ii+1]./λ)
            end
        end
    end

    #we will break into periodic and preperiodic cases for the last leg
    if S.orbit.preperiod == 0 #periodic
        u = sqrt(S.legs[1][1]./λ)
        thetau = (angle(u)/(2*pi)+1)%1
        
        if abs2(thetau - a/2) < abs2(thetau - (a+1)/2)
            if cregion == 'A'
                if S.kneading_sequence.items[end] == '2'
                    newLegs[end] = path_sqrt(S.legs[1]./λ)
                else
                    newLegs[end] = -1 .*path_sqrt(S.legs[1]./λ)
                end
            else
                if S.kneading_sequence.items[end] == '1'
                    newLegs[end] = path_sqrt(S.legs[1]./λ)
                else
                    newLegs[end] = -1 .*path_sqrt(S.legs[1]./λ)
                end
            end
        else
            if cregion == 'A'
                if S.kneading_sequence.items[end] == '1'
                    newLegs[end] = path_sqrt(S.legs[1]./λ)
                else
                    newLegs[end] = -1 .*path_sqrt(S.legs[1]./λ)
                end
            else
                if S.kneading_sequence.items[end] == '2'
                    newLegs[end] = path_sqrt(S.legs[1]./λ)
                else
                    newLegs[end] = -1 .*path_sqrt(S.legs[1]./λ)
                end
            end
        end

    else #preperiodic
        p = S.preperiod 
        u = sqrt(S.legs[p+1][1]./λ)
        thetau = (angle(u)/(2*pi)+1)%1

        if a/2 < thetau && thetau < (a+1)/2 #thetau is in the connected region
            if S.kneading_sequence.items[end] == cregion #then this is the correct preimage
                newLegs[end] = path_sqrt(S.legs[p+1]./λ)
            else
                newLegs[end] = -1 .* path_sqrt(S.legs[p+1]./λ)
            end
        else #thetau is in the disconnected region
            if S.kneading_sequence.items[end] == dregion #then this is the correct preimage
                newLegs[end] = path_sqrt(S.legs[p+1]./λ)
            else
                newLegs[end] = -1 .* path_sqrt(S.legs[p+1]./λ)
            end
        end
    end

    for leg in newLegs
        leg .+= (-1.0+0.0im)
    end
    newLegs .*= (2.0+0.0im)

    grow!(newLegs,10,10)
    for ii in eachindex(S.legs)
        S.legs[ii] = copy(newLegs[ii]) #done in this way so we can mutate S without mutating S
    end
    return S
end
        
function spideriterates(S0::Spider,n_iter::Int)
    list = [deepcopy(S0)]

    for i in 1:n_iter
        push!(list,deepcopy(mapspider!(S0)))
    end
    return list
end