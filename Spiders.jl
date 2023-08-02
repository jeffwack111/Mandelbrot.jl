using IterTools
include("SpiderFuncs.jl")
include("Sequences.jl")

struct SpiderLeg
    leg::Vector{ComplexF64}
    angle::Rational
end

function standard_leg(angle::Rational)
    r = collect(LinRange(1,100,1000))  #NOTE - it may be better to keep this as a 'linrange' but I don't understand what that means
    return SpiderLeg(r.*exp(1.0im*2*pi*angle),angle)
end

function grow!(S::SpiderLeg,scale::Real,num::Int)
    #Get the arrow pointing from the second to last point to the last point
    arrow = scale*(leg[end]-leg[end-1])
    #add new elements which are the last element plus the arrow
    for i in 1:num
        append!(S.leg,S.leg[end]+arrow)
    end
end

function prune!()
end

function info(SL::SpiderLeg)
    println("~~~~~~~~~~~~~~~~~~~~")
    λ = SL.leg[2][1]
    println(λ/2)
    points = length(SL.leg)*length(SL.leg[1])
    println(points)
    radius = sqrt(abs2(SL.leg[1][end]))
    println(radius)
end

#=
function spider_map(S::SpiderLegs,K::PString)

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