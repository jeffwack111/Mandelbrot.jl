using GLMakie
include("parameters/JuliaSet.jl")
include("interactive.jl")

fig = Figure()
ax1 = Axis(fig[1,1],aspect = 1)
ax2 = Axis(fig[1,2], aspect = 1)

colors = colorax(ax1)
juliaframe!(ax2,-1.0+0.0im,colors)
fig