using CairoMakie
include("CircleDynamics.jl")

fig = Figure()
ax = Axis(fig[1, 1])
r = 1.2
xlims!(-r,r)
ylims!(-r,r)

theta = 1//19    

orb = orbit(theta)

n = length(orb)

angle_list = orb[1:n+1]
x_list = cos.(2 .* pi .* angle_list)
y_list = sin.(2 .* pi .* angle_list)

lines!(x_list,y_list,linewidth = 2,color = "black")
fig


