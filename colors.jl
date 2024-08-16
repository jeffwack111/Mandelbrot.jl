using ColorSchemes
using Colors
using CairoMakie

n = 4
colors = [Colors.HSL{Float64}(360*i/n, 1, 0.5) for i in 0:(n-1)]
print(colors)
data = ones(n)

f, ax, plt = pie(data,
                 color = colors,
                 radius = 4,
                 inner_radius = 2,
                 strokecolor = :white,
                 strokewidth = 5,
                 axis = (autolimitaspect = 1, )
                )

f

