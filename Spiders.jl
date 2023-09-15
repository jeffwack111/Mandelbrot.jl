using IterTools
using CairoMakie
using ColorSchemes
include("SpiderFuncs.jl")
include("Sequences.jl")

struct SpiderInfo
    orbit::Vector{Rational}
    kneading_sequence::Vector{Char}
    preperiod::Int
end

function SpiderInfo(angle::Rational)

    orb = orbit(angle)
    K = kneadingsequence(orb)

    return SpiderInfo(orb.items,K.items,K.preperiod)

end

function standardlegs(SpIf::SpiderInfo)
    r = collect(LinRange(1,100,1000))  #NOTE - it may be better to keep this as a 'linrange' but I don't understand what that means
    legs  = Vector{ComplexF64}[]
    for theta in SpIf.orbit
        push!(legs,(cos(theta*2*pi)+1.0im*sin(theta*2*pi)) .* r)
    end
    legs[1] = legs[1] .- legs[1][1]
    return legs
end

function standardlegs(angle::Rational)
    info = SpiderInfo(angle)
    return standardlegs(info)
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

function stats(legs::Vector{Vector{ComplexF64}})
    println("~~~~~~~~~~~~~~~~~~~~")
    λ = legs[2][1]
    println(λ/2)
    points = length(leg)*length(leg[1])
    println(points)
    radius = sqrt(abs2(leg[1][end]))
    println(radius)
end

function devspidermap(angle::Rational)
    return devspidermap(SpiderInfo(angle),standardlegs(angle))
end


function devspidermap(S::SpiderInfo,Legs::Vector{Vector{ComplexF64}})

    fig = Figure()
    ax1 = Axis(fig[1, 1])
    r = 4
    xlims!(-r,r)
    ylims!(-r,r)

    n = length(Legs)

    for (j ,leg) in enumerate(Legs)
        lines!(real(leg),imag(leg),color = get(ColorSchemes.viridis, float(j)/float(n)))
    end

    λ = Legs[2][1]
    scaledLegs = Legs ./ λ

    ax2 = Axis(fig[1, 2])
    r = 4
    xlims!(-r,r)
    ylims!(-r,r)
    for (j ,leg) in enumerate(scaledLegs)
        lines!(real(leg),imag(leg),color = get(ColorSchemes.viridis, float(j)/float(n)))
    end

    ax3 = Axis(fig[2, 1])
    r = 4
    xlims!(-r,r)
    ylims!(-r,r)

    n = length(S.orbit)
    newLegs = Vector{Vector{ComplexF64}}(undef,n)

    regions = "AB"
    #First n-1 legs
    for ii in 1:n-1
        #This foot is guaranteed to be in the right half-plane thanks to the IEEE definition of sqrt
        #What region is it in? The right preimage of a foot is in region B iff the line connecting the foot to foot 2 crosses through leg 1
        line = (Legs[2][1],Legs[ii+1][1])
        pairs = partition(Legs[1],2,1)
        intersections = Vector{Bool}(undef,length(pairs))
        for (j,pair) in enumerate(pairs)
            intersections[j] = test_intersection(pair...,line...)
        end
        region_index = sum(intersections)%2 + 1
        if regions[region_index] == S.kneading_sequence[ii]
            newLegs[ii] = path_sqrt(Legs[ii+1]./λ)
        else
            newLegs[ii] = -1 .* path_sqrt(Legs[ii+1]./λ)
        end
    end

    line = (Legs[2][1],Legs[S.preperiod+1][1])
    pairs = partition(Legs[1],2,1)
    intersections = Vector{Bool}(undef,length(pairs))
    for (j,pair) in enumerate(pairs)
        intersections[j] = test_intersection(pair...,line...)
    end
    region_index = sum(intersections)%2 + 1
    if regions[region_index] == S.kneading_sequence[end]
        newLegs[end] = path_sqrt(Legs[S.preperiod+1]./λ)
    else
        newLegs[end] = -1 .* path_sqrt(Legs[S.preperiod+1]./λ)
    end

    for (j ,leg) in enumerate(newLegs)
        lines!(real(leg),imag(leg),color = get(ColorSchemes.viridis, float(j)/float(n)))
    end

    ax4 = Axis(fig[2, 2])
    r = 4
    xlims!(-r,r)
    ylims!(-r,r)

    for leg in newLegs
        leg .+= (-1.0+0.0im)
    end
    newLegs .*= (2.0+0.0im)
    grow!(newLegs,10,10)

    for (j ,leg) in enumerate(newLegs)
        lines!(real(leg),imag(leg),color = get(ColorSchemes.viridis, float(j)/float(n)))
    end

    return fig, newLegs

end

function spidermap(S::SpiderInfo,Legs::Vector{Vector{ComplexF64}})
    λ = Legs[2][1]

    regions = "AB"
    n = length(S.orbit)

    newLegs = Vector{Vector{ComplexF64}}(undef,n)

    #First n-1 legs
    for ii in 1:n-1
        #This foot is guaranteed to be in the right half-plane thanks to the IEEE definition of sqrt
        #What region is it in? The right preimage of a foot is in region B iff the line connecting the foot to foot 2 crosses through leg 1
        line = (Legs[2][1],Legs[ii+1][1])
        pairs = partition(Legs[1],2,1)
        intersections = Vector{Bool}(undef,length(pairs))
        for (j,pair) in enumerate(pairs)
            intersections[j] = test_intersection(pair...,line...)
        end
        region_index = sum(intersections)%2 + 1
        if regions[region_index] == S.kneading_sequence[ii]
            newLegs[ii] = path_sqrt(Legs[ii+1]./λ)
        else
            newLegs[ii] = -1 .* path_sqrt(Legs[ii+1]./λ)
        end
    end

    line = (Legs[2][1],Legs[S.preperiod+1][1])
    pairs = partition(Legs[1],2,1)
    intersections = Vector{Bool}(undef,length(pairs))
    for (j,pair) in enumerate(pairs)
        intersections[j] = test_intersection(pair...,line...)
    end
    region_index = sum(intersections)%2 + 1
    if regions[region_index] == S.kneading_sequence[end]
        newLegs[end] = path_sqrt(Legs[S.preperiod+1]./λ)
    else
        newLegs[end] = -1 .* path_sqrt(Legs[S.preperiod+1]./λ)
    end

    for leg in newLegs
        leg .+= (-1.0+0.0im)
    end
    newLegs .*= (2.0+0.0im)

    grow!(newLegs,10,10)

    return newLegs
    
end


