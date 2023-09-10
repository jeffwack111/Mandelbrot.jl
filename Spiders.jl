using IterTools
using CairoMakie
using ColorSchemes
include("SpiderFuncs.jl")
include("CircleDynamics.jl")

struct SpiderInfo
    orbit::Vector{Rational}
    kneading_sequence::String
    mapsto::Vector{Int}
end

function SpiderInfo(angle::Rational)
    
    a = angle/2
    b = (angle+1)/2
    self_itinerary = Char[]
    orb,mapsto = orbit(angle)

    for theta in orb
        if theta == a
            push!(self_itinerary,'a')
        elseif theta == b
            push!(self_itinerary,'b')
        elseif theta > a && theta < b
            push!(self_itinerary,'A')
        else
            push!(self_itinerary,'B')
        end
    end

    return SpiderInfo(orb,join(self_itinerary),mapsto)

end

function standard_legs(angle::Rational)
    info = SpiderInfo(angle)
    r = collect(LinRange(1,100,1000))  #NOTE - it may be better to keep this as a 'linrange' but I don't understand what that means
    legs  = Vector{ComplexF64}[]
    for theta in info.orbit
        push!(legs,(cos(theta*2*pi)+1.0im*sin(theta*2*pi)) .* r)
    end
    legs[1] = legs[1] .- legs[1][1]
    return legs
end

function grow!(legs::Vector{Vector{ComplexF64}},scale::Real,num::Int)
    #Get the arrow pointing from the second to last point to the last point
    for leg in legs
        arrow = scale*(leg[end]-leg[end-1])
        #add new elements which are the last element plus the arrow
        for i in 1:num
            append!(leg,leg[end]+arrow)
        end
    end
end


function prune!()
end

function spider_info(legs::Vector{Vector{ComplexF64}})
    println("~~~~~~~~~~~~~~~~~~~~")
    λ = legs[2][1]
    println(λ/2)
    points = length(leg)*length(leg[1])
    println(points)
    radius = sqrt(abs2(leg[1][end]))
    println(radius)
end

function dev_spider_map(angle::Rational)
    return dev_spider_map(SpiderInfo(angle),standard_legs(angle))
end


function dev_spider_map(S::SpiderInfo,Legs::Vector{Vector{ComplexF64}})

    fig = Figure()
    ax1 = Axis(fig[1, 1])
    r = 4
    xlims!(-r,r)
    ylims!(-r,r)

    newLegs = copy(Legs)

    n = length(newLegs)

    for (j ,leg) in enumerate(newLegs)
        lines!(real(leg),imag(leg),color = get(ColorSchemes.viridis, float(j)/float(n)))
    end


    ax2 = Axis(fig[1, 2])
    r = 4
    xlims!(-r,r)
    ylims!(-r,r)

    λ = Legs[2][1]
    newLegs = newLegs ./ λ

    for (j ,leg) in enumerate(newLegs)
        lines!(real(leg),imag(leg),color = get(ColorSchemes.viridis, float(j)/float(n)))
    end

    ax3 = Axis(fig[2, 1])
    r = 4
    xlims!(-r,r)
    ylims!(-r,r)
    newnewLegs = Vector{Vector{ComplexF64}}(undef,length(S.mapsto))

    regions = "AB"

    for (index,source) in enumerate(S.mapsto)
        #This foot is guaranteed to be in the right half-plane thanks to the IEEE definition of sqrt
        #What region is it in? The right preimage of a foot is in region B iff the line connecting the foot to foot 2 crosses through leg 1
        line = (newLegs[2][1],newLegs[source][1])
        pairs = partition(newLegs[1],2,1)
        intersections = Vector{Bool}(undef,length(pairs))
        for (j,pair) in enumerate(pairs)
            intersections[j] = test_intersection(pair...,line...)
        end
        region_index = sum(intersections)%2 + 1
        if regions[region_index] == S.kneading_sequence[index]
            newnewLegs[index] = path_sqrt(Legs[source]./λ)
        else
            newnewLegs[index] = -1 .* path_sqrt(Legs[source]./λ)
        end
        
    end
    
    for (j ,leg) in enumerate(newnewLegs)
        lines!(real(leg),imag(leg),color = get(ColorSchemes.viridis, float(j)/float(n)))
    end

    ax4 = Axis(fig[2, 2])
    r = 4
    xlims!(-r,r)
    ylims!(-r,r)

    for leg in newnewLegs
        leg .+= (-1.0+0.0im)
    end
    newnewLegs .*= (2.0+0.0im)

    for (j ,leg) in enumerate(newnewLegs)
        lines!(real(leg),imag(leg),color = get(ColorSchemes.viridis, float(j)/float(n)))
    end

    return fig

end

function spider_map(S::SpiderInfo,Legs::Vector{Vector{ComplexF64}})
    newLegs = Vector{ComplexF64}[]

    regions = "AB"

    for (index,source) in enumerate(S.mapsto)
        #This foot is guaranteed to be in the right half-plane thanks to the IEEE definition of sqrt
        #What region is it in? The right preimage of a foot is in region B iff the line connecting the foot to foot 2 crosses through leg 1
        λ = Legs[2][1]
        line = (λ,Legs[source][1])
        pairs = partition(L[1],2,1)
        intersections = Vector{Bool}(undef,length(pairs))
        for (j,pair) in enumerate(pairs)
            intersections[j] = test_intersection(pair...,line...)
        end
        region_index = sum(intersections)%2 + 1
        if regions[region_index] == S.kneading_sequence[index]
            push!(newLegs, 2.0 .* (path_sqrt(Legs[source]./λ) .- 1.0))
        else
            push!(newLegs, 2.0 .* (-1 .* path_sqrt(Legs[source]./λ) .- 1.0))
        end
        
    end

    grow!(newLegs,10,10)

    return newLegs
    
end


#=


    #first calculate the two lifts of the 1st leg

    first_leg = S.legs[1]
    boundary_one = lift_path(first_leg,λ,1)
    boundary_two = lift_path(first_leg,λ,-1)

    #We may need boundary_one later if this spider is periodic
    boundary = copy(boundary_one)

    #join them together, creating the boundary between the regions A and B
    popfirst!(boundary)
    append!(reverse!(boundary),boundary_two)


    n = length(K.itinerary)

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
            if K.itinerary[end] == '1'
                push!(newLegs,boundary_one) 
            else
                push!(newLegs,boundary_two)
            end
        else
            if K.itinerary[end] == '1'
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
=#