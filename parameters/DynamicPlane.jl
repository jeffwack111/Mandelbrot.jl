include("DynamicRays.jl")
include("RenderFractal.jl")
include("../spiders/Spiders.jl")
include("../trees/Lamination.jl")

function characteristicrays(theta::Rational) 

    c = spideriterate(theta,100)

    n = period(theta)
    orbit = [0.0+0.0im]
    for ii in 1:n-1
        push!(orbit,orbit[end]^2+c)
    end
    
    thetaprime = conjugate(theta)

    rays = dynamicrays(c,theta,100,10,20)
    conjrays = dynamicrays(c,thetaprime,100,10,20)

    julia = inverseiterate(c,22)

    fig = Figure()
    ax = Axis(fig[1, 1],aspect = 1)
    r = 2
    xlims!(-r,r)
    ylims!(-r,r)
    hidedecorations!(ax)
    #hidespines!(ax)

    plotrays!(ax,rays)
    plotrays!(ax,conjrays)

    scatter!(ax,real(julia),imag(julia),markersize = 1,color = "black")
    #scatter!(ax,real(orbit),imag(orbit))

    return fig
end