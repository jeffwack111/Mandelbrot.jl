using GLMakie
using Colors

n_colors = 20

activebox = Observable{Int}(1)

fig = Figure()
ax = Axis(fig[1,1],aspect = 1)
hidedecorations!(ax)
hidespines!(ax)
deactivate_interaction!(ax, :rectanglezoom)
deactivate_interaction!(ax, :dragpan)
deactivate_interaction!(ax, :scrollzoom)

sg = SliderGrid(fig[2,1],
                (label = "R", range = 0:0.01:1, startvalue = 1),
                (label = "G", range = 0:0.01:1, startvalue = 0),
                (label = "B", range = 0:0.01:1, startvalue = 0))

red = sg.sliders[1].value
green = sg.sliders[2].value
blue = sg.sliders[3].value

colorobservables = Observable([Colors.RGB(rand(3)...) for ii in 1:n_colors])
data = ones(n_colors)
inn = 0.4
out = 1

register_interaction!(ax, :select) do event::MouseEvent, ax
    if event.type === MouseEventTypes.leftdown
        ang = atan(event.data...)
        rad = sum([x*x for x in event.data])
        if rad < out*out && rad > inn*inn
            idx = ceil((-ang / pi + 1)*n_colors/2)
            activebox[] = idx
        end
    end
end

on(activebox) do activebox
    color = colorobservables[][activebox]
    set_close_to!(sg.sliders[1], Float64(color.r))
    set_close_to!(sg.sliders[2], Float64(color.g))
    set_close_to!(sg.sliders[3], Float64(color.b))
end

on(red) do red
    color = colorobservables[][activebox[]]
    colorobservables[][activebox[]] = Colors.RGB(red,color.g,color.b)
    notify(colorobservables)
end

on(green) do green
    color = colorobservables[][activebox[]]
    colorobservables[][activebox[]] = Colors.RGB(color.r,green,color.b)
    notify(colorobservables)
end

on(blue) do blue
    color = colorobservables[][activebox[]]
    colorobservables[][activebox[]] = Colors.RGB(color.r,color.g,blue)
    notify(colorobservables)
end

pie!(ax,data,color = colorobservables,radius = out,inner_radius = inn, offset = 3*π/2)
notify(activebox)

fig