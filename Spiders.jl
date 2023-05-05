import Base
using IterTools
include("SpiderFuncs.jl")

#Do we need PrePeriodicKneadingSequence and PeriodicKneadingSequence as separate types?
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

        new(orb, theta_itinerary, joint)
    end
end

function orbit(theta::Rational)
    orb = Rational{Int64}[]
    angle = theta
        while isempty(findall(x->x==angle,orb))
            push!(orb,angle)
            angle = angle*2
            angle = angle%1//1
        end
    joint = findall(x->x==angle,orb)[1]
    return orb,joint
end

function itinerary(theta::Real,orbit::Vector{<:Real})
    star_2 = theta/2
    star_1 = (theta+1)/2
    theta_itinerary = Any[]
    for angle in orbit
        if angle == star_1
            push!(theta_itinerary,1)
        elseif angle == star_2
            push!(theta_itinerary,2)
        elseif angle > star_2 && angle < star_1
            push!(theta_itinerary,"A")
        else
            push!(theta_itinerary,"B")
        end
    end
    return theta_itinerary
end

function standardlegs(orbit::Vector{<:Real})
    legs = Vector{ComplexF64}[]
    r = collect(LinRange(1,1000,1000))  #NOTE - it may be better to keep this as a 'linrange' but I don't understand what that means
    for angle in orbit
        push!(legs,r.*exp(1.0im*2*pi*angle).-exp(1.0im*2*pi*orbit[1]))
    end
    return legs
end

function spider_map(legs::Vector{Vector{ComplexF64}},K::KneadingSequence)

    newLegs = Vector{ComplexF64}[]

    λ = legs[2][1]

    #first calculate the two lifts of the 1st leg

    first_leg = legs[1]
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
        z = P_inv(legs[i+1][1],λ,1)
        #the function 'region' tells us which region the point is in
        #the kneading sequence tells us which region the point should be in
        if region(z,boundary) == K.itinerary[i]
            push!(newLegs,lift_path(legs[i+1],λ,1))
        else
            push!(newLegs,lift_path(legs[i+1],λ,-1))
        end
    end

    #Now we take care of the last leg, which depends on preperiodic or not
    if K.joint == 1
        #Periodic stuff
    else
        z = P_inv(legs[K.joint][1],λ,1)
        if region(z,boundary) == K.itinerary[n]
            push!(newLegs,lift_path(legs[K.joint],λ,1))
        else
            push!(newLegs,lift_path(legs[K.joint],λ,-1))
        end
    end

    return newLegs, boundary

end
