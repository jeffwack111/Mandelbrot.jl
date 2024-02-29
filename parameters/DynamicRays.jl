using CairoMakie
using ColorSchemes
include("../sequences/AngleDoubling.jl")
include("../spiders/SpiderFuncs.jl")

function dynamicrays(c::Complex,angle::Rational,R::Real,res::Int,depth::Int)
    orb = orbit(angle)

    raylist = Vector{ComplexF64}[]
    radii = collect(LinRange(R,sqrt(R),res))
    for theta in orb.items
        push!(raylist,exp(2im*pi*theta).*radii)
    end

    for jj in 1:depth
        for (target,source) in enumerate(goesto(orb))
            z = raylist[source][end-res+1]
            a = sqrt(-c + z)
            b = -sqrt(-c + z)

            da = abs2(raylist[target][end]-a)
            db = abs2(raylist[target][end]-b)

            if da < db
                append!(raylist[target],path_sqrt(-c .+ raylist[source][end-res+1:end]))
            elseif db < da
                append!(raylist[target],-1 .* path_sqrt(-c .+ raylist[source][end-res+1:end]))
            else
                return error("can't decide what to glue")
            end     
        end
    end
    return raylist
end

function landingpoint(c::Complex,angle::Rational,R::Real,res::Int,depth::Int)
    return dynamicrays(c::Complex,angle::Rational,R::Real,res::Int,depth::Int)[1][end]
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
    return plotrays(dynamicrays(c::Complex,angle::Rational,R::Real,res::Int,depth::Int))
end