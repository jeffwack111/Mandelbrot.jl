using GLMakie
using ColorSchemes

fig = Figure()
ax = Axis(fig[1, 1])
r = 4
xlims!(-r,r)
ylims!(-r,r)

n = length(list[end].legs)

record(fig,"test.gif",1:frames; framerate = 3) do i
    empty!(ax)
    spider = list[i]
    for (j ,leg) in enumerate(spider.legs)
        lines!(real(leg),imag(leg),color = get(ColorSchemes.viridis, float(j)/float(n)))
    end

    lines!(real(spider.boundary),imag(spider.boundary),linewidth = 4,color = "red")
end
