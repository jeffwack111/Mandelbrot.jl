using ColorSchemes
using Colors
using CairoMakie

n = 18
colors = [Colors.HSL{Float64}(30+360*i/n, 1, 0.5) for i in 0:(n-1)]
scene = Scene(size = (500,500))

function treecolors(A,B,distance)
    sat = 0.8
    light = 0.5
    return [Colors.HSL{Float64}(A, sat, light),
            Colors.HSL{Float64}(A+distance/4, sat, light),
            Colors.HSL{Float64}(A+distance/2, sat, light),
            Colors.HSL{Float64}(B-distance, sat, light),
            Colors.HSL{Float64}(B-distance/2, sat, light),
            Colors.HSL{Float64}(B, sat, light)]
end


function addsome(A,B,C)
    return [A,
            Colors.weighted_color_mean(0.125,C,A),
            Colors.weighted_color_mean(0.25,C,A),
            Colors.weighted_color_mean(0.25,C,B),
            Colors.weighted_color_mean(0.125,C,B),
            B]
end

cspace = XYZ
                
y = convert(cspace,Colors.RGB(1,1,0))
r = convert(cspace,Colors.RGB(1,0.2,0.2))
b = convert(cspace,Colors.RGB(0.2,0.2,1))


#colors = treecolors(0,220,70)

colors = addsome(r,b,y)

pie!(scene,ones(length(colors)),
                 color = colors,
                 radius = 1,
                 inner_radius = 0.1,
                 strokecolor = :white,
                 strokewidth = 0,
                )

scene


