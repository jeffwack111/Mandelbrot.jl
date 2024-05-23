using CairoMakie
using ColorSchemes
include("../sequences/AngleDoubling.jl")
include("../spiders/SpiderFuncs.jl")

## Conjecture: if the collection of rays is a sequence of rays, we will not need this function
function goesto(seq::Sequence)
    l = seq.preperiod
    k = period(seq)
    idx = append!(collect(1:l).+1,circshift(l+1:l+k,-1))
    return [(seq.items[x],seq.items[y]) for (x,y) in enumerate(idx)]
end 

### DANGER ###
#copy and pasted from HubbardTrees.jl
#similar to above function?
#=
function forwardimages(htree::Dict)
    return Dict([Pair(key,[shift(key)]) for key in keys(htree)])
end 
=#

function dynamicrays(c::Complex,angle::Rational,R::Real,res::Int,depth::Int)
    orb = orbit(angle)
    radii = collect(LinRange(R,sqrt(R),res))

    #similar to standard spider but is outer
    rays = Sequence{Vector{ComplexF64}}([exp(2im*pi*theta).*radii for theta in orb.items],orb.preperiod)

    for jj in 1:depth
        for (target,source) in goesto(rays)
            z = source[end-res+1]
            a = sqrt(-c + z)
            b = -sqrt(-c + z)

            da = abs2(target[end]-a)
            db = abs2(target[end]-b)

            if da < db
                append!(target,path_sqrt(-c .+ source[end-res+1:end]))
            elseif db < da
                append!(target,-1 .* path_sqrt(-c .+ source[end-res+1:end]))
            else
                return error("can't decide what to glue")
            end     
        end
    end
    return Dict([Pair(angle,ray) for (angle,ray) in zip(orb.items,rays.items)])
end


function landingpoint(c::Complex,angle::Rational,R::Real,res::Int,depth::Int)
    return dynamicrays(c::Complex,angle::Rational,R::Real,res::Int,depth::Int).items[1][end]
end

function plotrays!(ax,rays::Vector{Vector{ComplexF64}})
    n = length(rays)

    for (j ,ray) in enumerate(rays)
        lines!(ax,real(ray),imag(ray),color = get(ColorSchemes.rainbow, float(j)/float(n)))
    end

    return ax
end


function plotrays(rays::Vector{Vector{ComplexF64}})

    fig = Figure()
    ax = Axis(fig[1, 1],aspect = 1)
    r = 2
    xlims!(-r,r)
    ylims!(-r,r)

    n = length(rays)

    for (j ,ray) in enumerate(rays)
        lines!(ax,real(ray),imag(ray),color = get(ColorSchemes.cyclic_mygbm_30_95_c78_n256_s25, float(j)/float(n)))
    end

    return fig

end

function plotrays(c::Complex,angle::Rational,R::Real,res::Int,depth::Int)
    return plotrays(newdynamicrays(c::Complex,angle::Rational,R::Real,res::Int,depth::Int).items)
end