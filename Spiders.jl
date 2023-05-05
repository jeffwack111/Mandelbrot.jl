module Spiders

import Base
using IterTools

export Spider
export spider_map

abstract type Spider end

struct PeriodicSpider{T<:Complex} <: Spider
    angle::Rational
    legs::Vector{Vector{T}}
    itinerary::Vector{Any}

    function PeriodicSpider(angle::Rational,orbit::Vector{<:Real})
        legs = standardlegs(orbit)
        moveleg!(legs)

        itin = itinerary(angle,orbit)

        new{ComplexF64}(angle,legs,itin)
    end
end

struct PrePeriodicSpider{T<:Complex} <: Spider
    angle::Rational
    legs::Vector{Vector{T}}
    itinerary::Vector{Any}
    joint::Int

    function PrePeriodicSpider(angle::Rational,legs::Vector{Vector{ComplexF64}},itin::Vector{Any},joint::Int)
        new{ComplexF64}(angle,legs,itin,joint)
    end

    function PrePeriodicSpider(angle::Rational,orbit::Vector{<:Real},joint::Int)
        legs = standardlegs(orbit)
        moveleg!(legs)

        itin = itinerary(angle,orbit)

        new{ComplexF64}(angle,legs,itin,joint)
    end
end

function Spider(theta::Rational)
    orb, joint = orbit(theta)
    if joint == 1
        PeriodicSpider(theta,orb)
    else
        PrePeriodicSpider(theta,orb,joint)
    end
end



Base.show(io::IO, S::Spider) = display(hcat(S.itinerary, [round(leg[1],digits = 3) for leg in S.legs]))

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
    legs = Array{ComplexF64,1}[]
    r = collect(LinRange(1,50,100))  #NOTE - it may be better to keep this as a 'linrange' but I don't understand what that means
    for angle in orbit
        push!(legs,r.*exp(1.0im*2*pi*angle))
    end
    return legs
end

function moveleg!(spider)
    spider[1] = spider[1] .- spider[1][1]
end

function region(point::Complex,boundary::Vector{<:Complex})
    x2 = real(point)
    y2 = imag(point)
    
    n_segments = length(boundary)-1

    intersections = 0 

    for ii = 1:n_segments
        x3 = real(boundary[ii])
        y3 = imag(boundary[ii])
        
        x4 = real(boundary[ii+1])
        y4 = imag(boundary[ii+1])

        tn = y3*(x3-x4) - x3*(y3-y4)
        un = x3*y2-y3*x2
        d = y2*(y3-y4) - x2*(y3-y4)

        if sign(tn)*sign(d) == -1
            intersections += 0
        elseif sign(un)*sign(d) == -1
            intersections += 0
        elseif abs(un) == 0
            intersections +=0
        elseif abs(tn) > abs(d)
            intersections += 0
        elseif abs(un) > abs(d)
            intersections += 0
        else
            intersections += 1
        end
    end

    if intersections%2 == 0
        return "A"
    else
        return "B"
    end

end

function cross_cut(z1::Complex,z2::Complex)
    if sign(imag(z1)) == sign(imag(z2))
        return 1
    else
        if sign(real(z1)) == sign(real(z2))
            if real(z1) >= 0
                return 1
            else
                return -1
            end
        else
            #At this stage what we know is that both signs change
            #Now our points are sorted, so the line from a to b crosses the y-axis from left to right

            if real(z1) < real(z2)
                a = z1
                b = z2
            else
                a = z2
                b = z1
            end
            
            diff = b-a
            rise = imag(diff)
            run = real(diff)
            x = real(a)
            y = imag(b)

            #The line segment crosses the cut when the magnitude of the slope is larger than the slope of the line 
            #connecting the first point to the origin
            if abs(rise/run)>=abs(y/x)
                return -1
            else
                return 1
            end    
        end
    end
end

function lift_path(path::Vector{<:Complex},λ::Complex,branch)

    lift = [P_inv(path[1],λ,branch)]

    for pair in partition(path,2,1)
        #does the segment divided by lambda cross the cut?
        branch = branch * cross_cut(pair...)
        append!(lift,P_inv(pair[2],λ,branch))
    end
    
    return lift

end

function P_inv(z::Complex,λ::Complex,branch)
    return 2*(branch*sqrt(z/λ)-1)
end

function spider_map(spider::PrePeriodicSpider)

    newLegs = Vector{ComplexF64}[]

    λ = spider.legs[2][1]

    #first calculate the two lifts of the 1st leg

    first_leg = popfirst!(spider.legs)
    boundary = lift_path(first_leg,λ,1)
    boundary_two = lift_path(first_leg,λ,-1)

    #join them together, creating the boundary between the regions A and B
    popfirst!(boundary)
    append!(reverse!(boundary),boundary_two)

    #We now want to create an array of tuples
    #Each tuple will contain a leg to lift and the boundary to lift it into 

    push!(spider.legs,spider.legs[spider.joint-1])

    for (leg, reg) in zip(spider.legs, spider.itinerary)
        #first we will calculate an image of the foot, and determine if this point is in the correct region
        z = P_inv(leg[1],λ,1)
        #the function 'region' tells us which region the point is in
        #the kneading sequence tells us which region the point should be in
        if region(z,boundary) == reg
            push!(newLegs,lift_path(leg,λ,1))
        else
            push!(newLegs,lift_path(leg,λ,1))
        end
    end

    return PrePeriodicSpider(spider.angle,newLegs,spider.itinerary,spider.joint)

end

end