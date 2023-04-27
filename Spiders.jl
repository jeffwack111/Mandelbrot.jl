module Spiders

export Spider
export PeriodicSpider

abstract type Spider end

struct PeriodicSpider{T<:Complex} <: Spider
    angle::Rational
    legs::Vector{Vector{T}}
    itinterary::String
end

struct PrePeriodicSpider{T<:Complex} <: Spider
    angle::Rational
    legs::Vector{Vector{T}}
    itinterary::String
    joint::Int
end

function orbit(theta::Rational)
    itinerary = Rational{Int64}[]
    angle = theta
        while isempty(findall(x->x==angle,itinerary))
            push!(itinerary,angle)
            angle = angle*2
            angle = angle%1//1
        end
    period = length(itinerary) - findall(x->x==angle,itinerary)[1] + 1
    return itinerary,period
end

function Spider(theta::Rational)
