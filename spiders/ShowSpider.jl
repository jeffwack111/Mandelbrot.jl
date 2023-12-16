using CairoMakie
using ColorSchemes
include("Spiders.jl")

function spideriterates(angle::Rational,frames::Int)
    info = SpiderInfo(angle)

    list = [standardlegs(angle)]

    for i in 1:frames-1
        SL = list[end]
        push!(list,spidermap(info,SL))
    end
    return list
end

function showspider(angle::Rational,frames::Int)

    list = spideriterates(angle,frames)

    fig = Figure()
    ax = Axis(fig[1, 1])
    r = 8
    xlims!(-r,r)
    ylims!(-r,r)

    n = length(list[end])

    record(fig,"test.gif",1:frames; framerate = 3) do i
        empty!(ax)
        spider = list[i]
        for (j ,leg) in enumerate(spider)
            lines!(real(leg),imag(leg),color = get(ColorSchemes.viridis, float(j)/float(n)))
            text!(real(leg[1]),imag(leg[1]);text = "$j")
        end
        

    end 

    return list[end][2][end]/2
        
end 


