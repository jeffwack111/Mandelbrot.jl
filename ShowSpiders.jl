using GLMakie
using ColorSchemes

f = Figure()
ax = Axis(f[1, 1])
xlims!(-4,4)
ylims!(-4,4)

n = length(list[end].legs)

record(f,"test.gif",1:frames; framerate = 3) do i
    empty!(ax)
    spider = list[i]
    for (j ,leg) in enumerate(spider.legs)
        lines!(real(leg),imag(leg),color = get(ColorSchemes.viridis, float(j)/float(n)))
    end

    lines!(real(spider.boundary),imag(spider.boundary),linewidth = 4,color = "red")
end
