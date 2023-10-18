using CairoMakie
using ColorSchemes
include("AngleDoubling.jl")
include("SpiderFuncs.jl")

function dynamicrays(c::Complex,angle::Rational)
    orb = orbit(angle)

    raylist = Vector{ComplexF64}[]
    R = 100.0
    num = 100
    radii = collect(LinRange(R,sqrt(R),num))
    for theta in orb.items
        push!(raylist,exp(2im*pi*theta).*radii)
    end

    n = length(orb.items)

    for jj in 1:100
        for (target,source) in enumerate(goesto(orb))
            z = raylist[source][end-num+1]
            a = sqrt(-c + z)
            b = -sqrt(-c + z)

            da = abs2(raylist[target][end]-a)
            db = abs2(raylist[target][end]-b)

            if da < db
                append!(raylist[target],path_sqrt(-c .+ raylist[source][end-num+1:end]))
            elseif db < da
                append!(raylist[target],-1 .* path_sqrt(-c .+ raylist[source][end-num+1:end]))
            else
                return error("can't decide what to glue")
            end     
        end
    end
    return raylist
end

function plotrays(rays::Vector{Vector{ComplexF64}})

    fig = Figure()
    ax = Axis(fig[1, 1],aspect = 1)
    r = 2
    xlims!(-r,r)
    ylims!(-r,r)

    n = length(rays)

    for (j ,ray) in enumerate(rays)
        lines!(real(ray),imag(ray),color = get(ColorSchemes.viridis, float(j)/float(n)))
    end

    return fig

end

function plotrays(c::Complex,angle::Rational)
    return plotrays(dynamicrays(c::Complex,angle::Rational))
end