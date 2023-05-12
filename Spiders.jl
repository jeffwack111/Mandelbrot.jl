import Base
using IterTools
include("SpiderFuncs.jl")

#Shoulde we have PrePeriodicKneadingSequence and PeriodicKneadingSequence as separate types?
struct KneadingSequence
    orbit::Vector{Rational}
    itinerary::Vector{Any}
    joint::Int

    function KneadingSequence(theta::Rational)
        orb = Rational{Int64}[]

        star_2 = theta/2
        star_1 = (theta+1)/2
        theta_itinerary = Any[]

        angle = theta
        while isempty(findall(x->x==angle,orb))
            if angle == star_1
                push!(theta_itinerary,1)
            elseif angle == star_2
                push!(theta_itinerary,2)
            elseif angle > star_2 && angle < star_1
                push!(theta_itinerary,"A")
            else
                push!(theta_itinerary,"B")
            end

            push!(orb,angle)

            angle = angle*2
            angle = angle%1//1
        end

        joint = findall(x->x==angle,orb)[1]

        return new(orb, theta_itinerary, joint)
    end
end

struct SpiderLegs{T<:Complex}
    legs::Vector{Vector{T}}
    boundary::Vector{ComplexF64}

    #Standard constructor
    function SpiderLegs(legs::Vector{Vector{ComplexF64}},boundary::Vector{ComplexF64})
        return new{ComplexF64}(legs,boundary)
    end

    #Constructor which creates "default legs" from an orbit
    function SpiderLegs(orbit::Vector{<:Real})
        legs = Vector{ComplexF64}[]
        r = collect(LinRange(1,1000,1000))  #NOTE - it may be better to keep this as a 'linrange' but I don't understand what that means
        for angle in orbit
            push!(legs,r.*exp(1.0im*2*pi*angle).-exp(1.0im*2*pi*orbit[1]))
        end
        boundary::Vector{ComplexF64} = []
        return new{ComplexF64}(legs,boundary)
    end

end

function grow!(S::SpiderLegs,scale::Real,num::Int)
    for leg in S.legs
        #Get the arrow pointing from the second to last point to the last point
        arrow = scale*(leg[end]-leg[end-1])
        #add new elements which are the last element plus the arrow
        for i in 1:num
            append!(leg,leg[end]+arrow)
        end
    end
end

function prune!()
end

function info(SL::SpiderLegs)
    println("~~~~~~~~~~~~~~~~~~~~")
    λ = SL.legs[2][1]
    println(λ/2)
    len = length(SL.legs[1])
    println(len)
    radius = sqrt(abs2(SL.legs[1][end]))
    println(radius)
end

function spider_map(S::SpiderLegs,K::KneadingSequence)

    newLegs = Vector{ComplexF64}[]

    λ = S.legs[2][1]

    #first calculate the two lifts of the 1st leg

    first_leg = S.legs[1]
    boundary_one = lift_path(first_leg,λ,1)
    boundary_two = lift_path(first_leg,λ,-1)

    #We may need boundary_one later if this spider is periodic
    boundary = copy(boundary_one)

    #join them together, creating the boundary between the regions A and B
    popfirst!(boundary)
    append!(reverse!(boundary),boundary_two)


    n = length(K.orbit)

    for i in 1:n-1
        #first we will calculate an image of the foot, and determine if this point is in the correct region
        z = P_inv(S.legs[i+1][1],λ,1)
        #the function 'region' tells us which region the point is in
        #the kneading sequence tells us which region the point should be in
        if region(z,boundary) == K.itinerary[i]
            push!(newLegs,lift_path(S.legs[i+1],λ,1))
        else
            push!(newLegs,lift_path(S.legs[i+1],λ,-1))
        end
    end

    #Now we take care of the last leg, which depends on preperiodic or not
    if K.joint == 1
        theta_one = atan(real(boundary_one[end]),imag(boundary_one[end]))
        theta_two = atan(real(boundary_two[end]),imag(boundary_two[end]))
        if theta_one > theta_two
            if K.itinerary[end] == 1
                push!(newLegs,boundary_one) 
            else
                push!(newLegs,boundary_two)
            end
        else
            if K.itinerary[end] == 1
                push!(newLegs,boundary_two) 
            else
                push!(newLegs,boundary_one)
            end
        end  
    else
        z = P_inv(S.legs[K.joint][1],λ,1)
        if region(z,boundary) == K.itinerary[n]
            push!(newLegs,lift_path(S.legs[K.joint],λ,1))
        else
            push!(newLegs,lift_path(S.legs[K.joint],λ,-1))
        end
    end

    SL = SpiderLegs(newLegs,boundary)

    grow!(SL,10,10)
    info(SL)

    return SL

end