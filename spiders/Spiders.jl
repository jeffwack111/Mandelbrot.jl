using IterTools
include("SpiderFuncs.jl")
include("../sequences/AngleDoubling.jl")

struct SpiderInfo
    orbit::Vector{Rational}
    kneading_sequence::Vector{Char}
    preperiod::Int
end
    

function SpiderInfo(theta::Rational)

    orb = orbit(theta)
    K = kneadingsequence(theta)

    return SpiderInfo(orb.items,K.items,K.preperiod)

end

function standardlegs(SpIf::SpiderInfo)
    r = collect(LinRange(100,1,10))  #NOTE - it may be better to keep this as a 'linrange' but I don't understand what that means
    legs  = Vector{ComplexF64}[]
    for theta in SpIf.orbit
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

function spidermap(S::SpiderInfo, Legs::Vector{Vector{ComplexF64}})
    λ = Legs[2][end] #the parameter of our polynomial map z -> λ(1+z/2)^2
    a = (angle(Legs[1][1]/λ)/(2*pi)+1)%1 #the angle whose halves will define the boundary between regions A and B at infinity. in full turns

    n = length(S.orbit)

    newLegs = Vector{Vector{ComplexF64}}(undef,n)

    #leg 1 goes first
    newLegs[1] = path_sqrt(Legs[2]./λ)

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
        u = sqrt(Legs[ii+1][1]./λ)
        thetau = (angle(u)/(2*pi)+1)%1

        if a/2 < thetau && thetau < (a+1)/2 #thetau is in the connected region
            if S.kneading_sequence[ii] == cregion #then this is the correct preimage
                newLegs[ii] = path_sqrt(Legs[ii+1]./λ)
            else
                newLegs[ii] = -1 .* path_sqrt(Legs[ii+1]./λ)
            end
        else #thetau is in the disconnected region
            if S.kneading_sequence[ii] == dregion #then this is the correct preimage
                newLegs[ii] = path_sqrt(Legs[ii+1]./λ)
            else
                newLegs[ii] =  -1 .* path_sqrt(Legs[ii+1]./λ)
            end
        end
    end

    #we will break into periodic and preperiodic cases for the last leg
    if S.preperiod == 0 #periodic
        u = sqrt(Legs[1][1]./λ)
        thetau = (angle(u)/(2*pi)+1)%1
        
        if abs2(thetau - a/2) < abs2(thetau - (a+1)/2)
            if cregion == 'A'
                if S.kneading_sequence[end] == '2'
                    newLegs[end] = path_sqrt(Legs[1]./λ)
                else
                    newLegs[end] = -1 .*path_sqrt(Legs[1]./λ)
                end
            else
                if S.kneading_sequence[end] == '1'
                    newLegs[end] = path_sqrt(Legs[1]./λ)
                else
                    newLegs[end] = -1 .*path_sqrt(Legs[1]./λ)
                end
            end
        else
            if cregion == 'A'
                if S.kneading_sequence[end] == '1'
                    newLegs[end] = path_sqrt(Legs[1]./λ)
                else
                    newLegs[end] = -1 .*path_sqrt(Legs[1]./λ)
                end
            else
                if S.kneading_sequence[end] == '2'
                    newLegs[end] = path_sqrt(Legs[1]./λ)
                else
                    newLegs[end] = -1 .*path_sqrt(Legs[1]./λ)
                end
            end
        end

    else #preperiodic
        p = S.preperiod 
        u = sqrt(Legs[p+1][1]./λ)
        thetau = (angle(u)/(2*pi)+1)%1

        if a/2 < thetau && thetau < (a+1)/2 #thetau is in the connected region
            if S.kneading_sequence[end] == cregion #then this is the correct preimage
                newLegs[end] = path_sqrt(Legs[p+1]./λ)
            else
                newLegs[end] = -1 .* path_sqrt(Legs[p+1]./λ)
            end
        else #thetau is in the disconnected region
            if S.kneading_sequence[end] == dregion #then this is the correct preimage
                newLegs[end] = path_sqrt(Legs[p+1]./λ)
            else
                newLegs[end] = -1 .* path_sqrt(Legs[p+1]./λ)
            end
        end
    end

    for leg in newLegs
        leg .+= (-1.0+0.0im)
    end
    newLegs .*= (2.0+0.0im)

    grow!(newLegs,10,10)

    return newLegs
end
        
function spideriterate(angle::Rational,frames::Int)
    info = SpiderInfo(angle)

    list = [standardlegs(angle)]

    for i in 1:frames-1
        SL = list[end]
        push!(list,spidermap(info,SL))
    end
    return list[end][2][end]/2
end