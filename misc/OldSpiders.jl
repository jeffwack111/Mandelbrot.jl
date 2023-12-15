#NOTE: These spider legs are indexed with their feet (cloesest to zero) in position 1, and their shoulders (closest to infinity) at the end
#This may be cause for calls to reverse()

function devHSspidermap(S::SpiderInfo,Legs::Vector{Vector{ComplexF64}})

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

function HSspidermap(S::SpiderInfo,Legs::Vector{Vector{ComplexF64}})
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