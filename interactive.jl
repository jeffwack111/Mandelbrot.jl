using GLMakie
using Colors

fig = Figure()

sg = SliderGrid(fig[2,1],
                (label = "R", range = 0:0.01:1, startvalue = 1),
                (label = "G", range = 0:0.01:1, startvalue = 0),
                (label = "B", range = 0:0.01:1, startvalue = 0))

sliderobservables = [s.value for s in sg.sliders]
col = lift(sliderobservables...) do r,g,b
    Colors.RGB(r,g,b)
end

Box(fig[1, 1], color = col, cornerradius = 20, strokevisible = false)

fig