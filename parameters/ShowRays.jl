using GLMakie
using ColorSchemes
include("DynamicRays.jl")

@recipe(RayPlot,ExtRays) do scene 
    Theme()
end

function Makie.plot!(myplot::RayPlot)
    rays = collect(values(myplot.ExtRays[].rays))
    n = length(rays)
    for (j ,ray) in enumerate(rays)
        lines!(myplot,real(ray),imag(ray),color = get(ColorSchemes.rainbow, float(j)/float(n)))
    end
    return myplot
end