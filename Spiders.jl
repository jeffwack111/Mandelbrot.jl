module Spiders

export Spider

abstract type Spider end

struct PeriodicSpider{T<:Complex} <: Spider
    angle::Rational
    legs::Vector{Vector{T}}
    itinterary::Vector{Any}

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
    itinterary::Vector{Any}
    joint::Int

    function PrePeriodicSpider(angle::Rational,orbit::Vector{<:Real},joint::Int)
        legs = standardlegs(orbit)
        moveleg!(legs)

        itin = itinerary(angle,orbit)

        new{ComplexF64}(angle,legs,itin,joint)
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

function Spider(theta::Rational)
    orb, joint = orbit(theta)
    if joint == 1
        PeriodicSpider(theta,orb)
    else
        PrePeriodicSpider(theta,orb,joint)
    end
end

end