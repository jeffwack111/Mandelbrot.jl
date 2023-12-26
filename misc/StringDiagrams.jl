using CairoMakie
include("..//sequences//AngleDoubling.jl")

denominator_list = [19, 29, 37, 53, 59, 61, 67, 83, 101, 107, 131, 139, 149, 163, 173, 179, 181, 197, 211, 227, 269, 293, 317, 347, 349, 373]

fig = Figure(backgroundcolor = (:tomato,0),resolution = (1000, 1000))
ax = Axis(fig[1, 1],aspect = 1,backgroundcolor = (:tomato,0))
r = 1.5
xlims!(-r,r)
ylims!(-r,r)
hidedecorations!(ax)
hidespines!(ax)

for denominator in denominator_list
    theta = 1//denominator

    orb = orbit(theta)
    n = length(orb.items)
    angle_list = orb[1:n+1]
    x_list = cos.(2 .* pi .* angle_list)
    y_list = sin.(2 .* pi .* angle_list)

    empty!(ax)
    lines!(x_list,y_list,linewidth = 1,color = "black")
    save("$denominator.png", fig)
end

